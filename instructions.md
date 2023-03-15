How to use the 3D Hough Transformations for Identifying Minimum Ionizing Particle Tracks
=========================================================================================
Part 1: Cloning the directory from GitHub

1) mkdir mip_tracking
2) cd mip_tracking
3) git clone https://github.com/pbutti/hough3d.git
4) cd hough3d
5) git checkout david ("david" branch has the changes that are implemented to the 3D Hough Transformations) 



Part 2: Installing Eigen and using MakeFile
1) go to the website (https://eigen.tuxfamily.org/index.php?title=Main_Page) and download Eigen 3.4.0 (make sure you know where the location of it is)
2) cd mip_tracking
3) mkdir EIGEN
4) cd EIGEN
5) move the eigen-3.4.0 directory into here: mv [full path location to eigen-3.4.0] .
6) cd ../hough3d
7) we need to make a small change into the MakeFile: change the line "LIBEIGEN = something" to "LIBEIGEN = [full path location to eigen-3.4.0]"
8) save the changes and then in the hough3d directory, execute: make

Part 3a: Using it on a sample 
* hough3d_v2.cpp is the script that contains the quality cuts necessary to identify MIP tracks. 
* hough3d.cpp is the default 3d line fitting script that will not be that useful.
* samples and root files need to be "reformatted" in order for them to be operated on by the hough transformations
-> this is because we only save a select number of variables and also create new variables that are not present in the existing LDMX files

1) cd mip_tracking
2) mkdir PROCESSING
3) cd PROCESSING
4) git clone git clone https://github.com/davidgjiang/MIP_Tracking.git
5) cd MIP_Tracking
6) any samples that you want to use needs to pass through skimmer.py 
-> make sure to edit the script accordingly, for example change the input file locations, filenames, the branch suffix from 'v3_v12' to whatever, etc.
-> just need to create some directory and put all of your samples in there. put that location into skimmer.py, then simply execute skimmer.py in the MIP_Tracking directory

Part 3b: Using it on a sample
1) cd mip_tracking/hough3d
2) open up the hough3d_v2.cpp script and look for the line that requires the "input file" and "output file"
3) make sure the input file location is correct and the output file location is what you want it to be
4) feel free to adjust the parameters and values of the cuts (should be line 190) or look for: "if (fracLayerGaps >= 0.4 || fracLayersMultHits >= 0.3 || numLayersWithHits <= 3 || maxInARow <= 2) {"
5) now to compile.

-> execute: g++ -I[location of root include] -L[location of root lib] -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint 
-lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -Wl,-rpath, [location of root lib] -stdlib=libc++ -lpthread -lm -ldl -std=c++14 -o 
./hough3d_v2.exe vector3d.o pointcloud.o hough.o sphere.o hough3d_v2.cpp

-> example (mine): g++ -I/usr/local/Cellar/root/6.26.06_1/include/root/ -L/usr/local/Cellar/root/6.26.06_1/lib/root -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d 
-lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -Wl,-rpath,/usr/local/Cellar/root/6.26.06_1/lib/root 
-stdlib=libc++ -lpthread -lm -ldl -std=c++14 -o ./hough3d_v2.exe vector3d.o pointcloud.o hough.o sphere.o hough3d_v2.cpp

6) Final step is to execute: ./hough3d_v2.exe -dx [some value] -minvotes [another value] -nlines [and another value]
-> example: ./hough3d_v2.exe -dx 10 -minvotes 3 -nlines 2

Now your output file will have the necessary information about the center point of each MIP Track line and its 2 direction vectors!





