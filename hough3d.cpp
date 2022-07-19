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

struct compare
{
	int key;
	compare(int const &i): key(i) { }

	bool operator()(int const &i) {
		return (i == key);
	}
};

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

// minimum distance between two 3D lines (each line defined by 2 points)
double dist3Dlines(Vector3d* start1,Vector3d* end1,Vector3d* start2,Vector3d* end2) {  
  
  // direction vectors
  Vector3d e1 = start1 - end1;
  Vector3d e2 = start2 - end2;

  // vector perpendicular to both lines
  Vector3d n;
  n.x = (e1.y)*(e2.z) - (e1.z)*(e2.y);
  n.y = (e1.z)*(e2.x) - (e1.x)*(e2.z);
  n.z = (e1.x)*(e2.y) - (e1.y)*(e2.x);

  // calculate distance
  Vector3d r = start1 - start2;
  double n_mag = sqrt((pow(n.x,2) + pow(n.y,2)+ pow(n.z,2)),0.5)
  
  return (n.x*r.x + n.y*r.y + n.z+r.z)/n_mag;
}

int main(int argc, char ** argv) {
  
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
  
  // output file
  TFile f1("skimmed_kaons_hough.root", "recreate"); 
  TTree t1("Events", "Events");

  // input file
  TFile *f = TFile::Open("skimmed_kaons/skimmed_kaons_0.root");
  f->cd();
  
  // get the tree from the file and assign it to a new local variable
  TTree *tree = (TTree*)(f->Get("Events")); 

  // declare local branch variables
  vector<double> * hitX = 0;
  vector<double> * hitY = 0;
  vector<double> * hitZ = 0;
  vector<int> * pdgID = 0;
  vector<double> * photonX = 0;
  vector<double> * photonY = 0;
  vector<double> * photonZ = 0;
  vector<double> * electronX = 0;
  vector<double> * electronY = 0;
  vector<double> * electronZ = 0;
  
  int nlines=0;

  // set addresses for tree branches
  tree->SetBranchAddress("hitX", &hitX);
  tree->SetBranchAddress("hitY", &hitY);
  tree->SetBranchAddress("hitZ", &hitZ);
  tree->SetBranchAddress("pdgID", &pdgID); 
  tree->SetBranchAddress("photonX", &photonX);
  tree->SetBranchAddress("photonY", &photonY);
  tree->SetBranchAddress("photonZ", &photonZ);
  tree->SetBranchAddress("electronX", &electronX);
  tree->SetBranchAddress("electronY", &electronY);
  tree->SetBranchAddress("electronZ", &electronZ);


  //Fit results
  std::vector<double> a_x; 
  std::vector<double> a_y;
  std::vector<double> a_z;
  std::vector<double> b_x;   
  std::vector<double> b_y;
  std::vector<double> b_z;
  std::vector<std::vector<double>> pdg_assoc;
  std::vector<double> distFromPTraj;
  std::vector<double> distFromETraj;
  

  // set the addresses for output tree
  //t1->Branch("hitN",&hitN,"/I");
  t1.Branch("hitX",&hitX);
  t1.Branch("hitY",&hitY);
  t1.Branch("hitZ",&hitZ);
  t1.Branch("pdgID",&pdgID);
  t1.Branch("photonX",&photonX);
  t1.Branch("photonY",&photonY);
  t1.Branch("photonZ",&photonZ);
  t1.Branch("electronX",&electronX);
  t1.Branch("electronY",&electronY);
  t1.Branch("electronZ",&electronZ);
  t1.Branch("nlines",&nlines,"/I");
  t1.Branch("ax",&a_x);
  t1.Branch("ay",&a_y);
  t1.Branch("az",&a_z);
  t1.Branch("bx",&b_x);
  t1.Branch("by",&b_y);
  t1.Branch("bz",&b_z);
  t1.Branch("pdg_assoc",&pdg_assoc);
  t1.Branch("distFromPTraj",&distFromPTraj);
  t1.Branch("distFromETraj",&distFromETraj);
  
  // loop over all events
  for (int event = 0; event < tree->GetEntries(); event++) {
    
    tree->GetEntry(event); // this fills our local variables that we created with the ith event

    Vector3d minP, maxP, minPshifted, maxPshifted; // bounding box of point cloud
    
    double d; // diagonal length of point cloud
    PointCloud X; // point cloud
    
    // loop over all hits and add them to the point cloud
    for (int hit = 0; hit < hitX->size(); hit++) {
      X.addPoint(X,hitX->at(hit), hitY->at(hit), hitZ->at(hit));
    }

    Vector3d pTrajStart, pTrajEnd, eTrajStart, eTrajEnd; // photon and electron starting and ending points

    for (int layer = 0; layer < photonX.size(); layer++) {
      // starting point of photon and electron trajectories
      if (layer == 0) {
        pTrajStart.x = photonX->at(layer);
        pTrajStart.y = photonY->at(layer);
        pTrajStart.z = photonZ->at(layer);
        eTrajStart.x = electronX->at(layer);
        eTrajStart.y = electronY->at(layer);
        eTrajStart.z = electronZ->at(layer);
      }
      // ending point of photon and electron trajectories
      if (layer == photonX.size() - 1) {
        pTrajEnd.x = photonX->at(layer);
        pTrajEnd.y = photonY->at(layer);
        pTrajEnd.z = photonZ->at(layer);
        eTrajEnd.x = electronX->at(layer);
        eTrajEnd.y = electronY->at(layer);
        eTrajEnd.z = electronZ->at(layer);
      }
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
    
    double sameLayerCounter = 0.0; // counter for number of layers with same layer hits 
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
      //std::cout<<"a = " << a << ", b = " << b << std::endl;

      // loop over all of the points in the point cloud Y (hits in the track)
      std::vector<double> pdgsOfTrack;
      std::vector<unsigned int> reference;
      for (unsigned int i=0; i<Y.points.size(); i++) {
        Vector3d p = Y.points[i] + X.shift;
        reference.push_back(i); // vector for remembering points that were used 

        // 1. loop over all the other points in the track (PointCloud Y)
        for (unsigned int j=0; j<Y.points.size(); j++) {

          // make sure to not double count points
          if (j != i && !any_of(reference.begin(),reference.end(),compare(j))) {
            Vector3d p2 = Y.points[j] + X.shift;
            // increment the counter if 2 points appear in the same layer 
            if (p.z == p2.z) {
              sameLayerCounter += 1;
            }
          }
        }

        pdgsOfTrack.clear(); // clear the pdg IDs associated with the track

        // 2. loop over all of the hits in the event to save pdg ID (only if track doesnt have same layer hits)
        for (unsigned int hit=0; hit<hitX->size(); hit++) {
          if (abs(hitX->at(hit) - p.x) < 1e-6 && abs(hitY->at(hit) - p.y) < 1e-6 && abs(hitZ->at(hit) - p.z) < 1e-6) {
            double id = (double) pdgID->at(hit);
            //cout << id << endl;
            pdgsOfTrack.push_back(id);
            break;
          }
        }
      }

      // save distance from fitted track to photon trajectory and electron trajectory
      Vector3d startPoint = a - b*100;
      Vector3d endPoint = a + b*100;
      double pDist = dist3Dlines(pTrajStart,pTrajEnd,startPoint,endPoint);
      double eDist = dist3Dlines(eTrajStart,eTrajEnd,startPoint,endPoint);

      pDistances.push_back(pDist);
      eDistances.push_back(eDist)

      if (pdgsOfTrack.size() != 0) {
        pdg_assoc.push_back(pdgsOfTrack);
      }
      pdgsOfTrack.clear();

      X.removePoints(Y);

      a_x.push_back(a.x);
      a_y.push_back(a.y);
      a_z.push_back(a.z);

      b_x.push_back(b.x); 
      b_y.push_back(b.y);
      b_z.push_back(b.z);

      distFromETraj.push_back(eDist);
      distFromPTraj.push_back(pDist);

      nlines++; // only increment the number of tracks if no hits in the same layer

      sameLayerCounter = 0;

    } while ((X.points.size() > 1) && 
            ((opt_nlines == 0) || (opt_nlines > nlines)));

    std::cout<<"Number of lines: "<<nlines<<std::endl;

    t1.Fill();

    nlines = 0; 
    a_x.clear();
    a_y.clear();
    a_z.clear();

    b_x.clear();
    b_y.clear();
    b_z.clear();
    pdg_assoc.clear();

    distFromETraj.clear();
    distFromPTraj.clear();

    // clean up
    delete hough;

  } // event loop

  f1.cd();
  t1.Write();
  f1.Close();

  return 0;
}