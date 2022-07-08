#include<cstdio>
#include<cstdlib>
#include<iostream>

#include "TTree.h"
#include "TFile.h"

#include "vector3d.h"
#include "pointcloud.h"
#include "hough.h"

#include <stdio.h>
#include <string.h>
#include <Eigen/Dense>

using Eigen::MatrixXf;
using namespace std;

// orthogonal least squares fit with libeigen
double orthogonal_LSQ(const PointCloud &pc, Vector3d* a, Vector3d* b){
  double rc = 0.0; // rc = largest eigenvalue

  // anchor point is mean value
  *a = pc.meanValue();

  // copy points to libeigen matrix
  Eigen::MatrixXf points = Eigen::MatrixXf::Constant(pc.points.size(), 3, 0);
  for (int i = 0; i < points.rows(); i++) {
    points(i,0) = pc.points.at(i).x;
    points(i,1) = pc.points.at(i).y;
    points(i,2) = pc.points.at(i).z;
  }

  // compute scatter matrix ...
  MatrixXf centered = points.rowwise() - points.colwise().mean();
  MatrixXf scatter = (centered.adjoint() * centered);

  // ... and its eigenvalues and eigenvectors
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eig(scatter);
  Eigen::MatrixXf eigvecs = eig.eigenvectors();

  // we need eigenvector to largest eigenvalue
  // libeigen yields it as LAST column
  b->x = eigvecs(0,2); b->y = eigvecs(1,2); b->z = eigvecs(2,2);
  rc = eig.eigenvalues()(2);

  return (rc);
}

int main(int argc, char ** argv) {
  // output file
  TFile f1("skimmed_events_pdg_hough.root", "recreate");
  TTree t1("Events", "Events");

  // input file
  TFile *f = TFile::Open("skimmed_events_pdg1.root");
  f->cd();
  
  // get the tree from the file and assign it to a new local variable
  TTree *tree = (TTree*)(f->Get("Events")); 

  // declare local branch variables
  vector<double> * hitX = 0;
  vector<double> * hitY = 0;
  vector<double> * hitZ = 0;
  vector<int> * pdgID = 0;
  int hitN;

  // set addresses for tree branches
  tree->SetBranchAddress("hitN", &hitN);
  tree->SetBranchAddress("hitX", &hitX);
  tree->SetBranchAddress("hitY", &hitY);
  tree->SetBranchAddress("hitZ", &hitZ);
  tree->SetBranchAddress("pdgID", &pdgID); 
  
  cout << "number of events: " << tree->GetEntries() << std::endl;
  
  // loop over all events
  for (int event = 0; event < tree->GetEntries(); event++) {
    
    // default values for command line options
    double opt_dx = 0.0;
    int opt_nlines = 0;
    int opt_minvotes = 0;

    // number of icosahedron subdivisions for direction discretization
    int granularity = 4;
    int num_directions[7] = {12, 21, 81, 321, 1281, 5121, 20481};

    // parse command line
    for (int i=1; i<argc; i++) {
      if (0 == strcmp(argv[i], "-dx")) {
        i++;
        if (i<argc) opt_dx = atof(argv[i]);
      }
      else if (0 == strcmp(argv[i], "-nlines")) {
        i++;
        if (i<argc) opt_nlines = atoi(argv[i]);
      }
      else if (0 == strcmp(argv[i], "-minvotes")) {
        i++;
        if (i<argc) opt_minvotes = atoi(argv[i]);
      }
    }

  

    tree->GetEntry(event); // this fills our local variables that we created with the ith event

    Vector3d minP, maxP, minPshifted, maxPshifted; // bounding box of point cloud
    
    double d; // diagonal length of point cloud
    PointCloud X; // point cloud
    
    // loop over all hits and add them to the point cloud
    for (int hit = 0; hit < hitX->size(); hit++) {
      X.addPoint(X,hitX->at(hit), hitY->at(hit), hitZ->at(hit));
    }

    // cannot make point clouds with < 2 points
    if (X.points.size() < 2) {
      fprintf(stderr, "Error: point cloud has less than two points\n");
      return 1;
    }

    // center the point cloud and compute a new bounding box
    X.getMinMax3D(&minP, &maxP);
    d = (maxP-minP).norm();
    if (d == 0.0) {
      fprintf(stderr, "Error: all points in point cloud identical\n");
      return 1;
    }

    X.shiftToOrigin();
    X.getMinMax3D(&minPshifted, &maxPshifted);

    // estimate size of Hough space
    if (opt_dx == 0.0) {
      opt_dx = d / 64.0;
    }
    else if (opt_dx >= d) {
      fprintf(stderr, "Error: dx too large\n");
      return 1;
    }
    double num_x = floor(d / opt_dx + 0.5);
    double num_cells = num_x * num_x * num_directions[granularity];

    // first Hough transform
    Hough* hough;
    try {
      hough = new Hough(minPshifted, maxPshifted, opt_dx, granularity);
    } catch (const std::exception &e) {
      fprintf(stderr, "Error: cannot allocate memory for %.0f Hough cells"
              " (%.2f MB)\n", num_cells, 
              (double(num_cells) / 1000000.0) * sizeof(unsigned int));
      return 2;
    }
    hough->add(X);

    // iterative Hough transform
    // (Algorithm 1 in IPOL paper)
    PointCloud Y;	// points close to line
    double rc;
    unsigned int nvotes;
    int nlines = 0;
    
    do {
      Vector3d a; // anchor point of line
      Vector3d b; // direction of line

      hough->subtract(Y); // do it here to save one call
      
      nvotes = hough->getLine(&a, &b);

      X.pointsCloseToLine(a, b, opt_dx, &Y);

      rc = orthogonal_LSQ(Y, &a, &b);
      if (rc==0.0) break;

      X.pointsCloseToLine(a, b, opt_dx, &Y);
      nvotes = Y.points.size();
      if (nvotes < (unsigned int)opt_minvotes) break;

      rc = orthogonal_LSQ(Y, &a, &b);
      if (rc==0.0) break;

      a = a + X.shift;
      std::cout<<a<<" " << b << std::endl;
      nlines++;

    } while ((X.points.size() > 1) && 
            ((opt_nlines == 0) || (opt_nlines > nlines)));
    

    std::cout<<"Number of lines:: "<<nlines<<std::endl;


    // clean up
    delete hough;
    
  }

  return 0;
}