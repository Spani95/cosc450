#include "SimplePly.h"
#include <iostream>
#include <stdlib.h>
#include <utility>
#include <math.h>
#include <map>

using namespace std;

double distance(vector<PlyPoint> plane, PlyPoint p);
Eigen::Vector3d crossProduct(Eigen::Vector3d p1, Eigen::Vector3d p2);
PlyPoint find_point(SimplePly ply, Eigen::Vector3i null_colour);
vector<pair<size_t, PlyPoint> > remove_points(vector<pair<size_t, PlyPoint> > plyPoints, size_t points);

int main (int argc, char *argv[]) {

  // Check the commandline arguments.
  if (argc != 5) {
    cout << "Usage: planeFinder <input file> <output file> <number of planes> <point-plane threshold> <number of RANSAC trials>" << std::endl;
    return -1;
  }
  int nPlanes = atoi(argv[3]);
  double threshold = atof(argv[4]);

  cout << "Searching for " << nPlanes << " planes" << endl;
  cout << "Using a point-plane threshold of " << threshold << " units" << endl;
  cout << "Applying RANSAC with k trials" << endl;  

  // Storage for the point cloud.
  SimplePly ply;

  // Read in the data from a PLY file
  cout << "Reading PLY data from " << argv[1] << endl;
  if (!ply.read(argv[1])) {
    cout << "Could not read PLY data from file " << argv[1] << endl;
    return -1;
  }
  cout << "Read " << ply.size() << " points" << endl;

  // Recolour points - here we are just doing colour based on index
  cout << "Recolouring points" << endl;
  vector<Eigen::Vector3i> colours;
  colours.push_back(Eigen::Vector3i(255,0,0));
  colours.push_back(Eigen::Vector3i(0,255,0));
  colours.push_back(Eigen::Vector3i(0,0,255));
  colours.push_back(Eigen::Vector3i(255,0,255));
  colours.push_back(Eigen::Vector3i(255,255,0));
  colours.push_back(Eigen::Vector3i(0,255,255));
  // Can add more colours as needed

  map<pair<unsigned long, PlyPoint> > plyPoints_copy;
  /*
  // Setting all points to a null colour (black)
  Eigen::Vector3i null_colour = Eigen::Vector3i(0,0,0);
  for (size_t ix = 0; ix < ply.size(); ++ix) {
    ply[ix].colour = null_colour;
  }
  */
  cout << "making copy" << endl;
  for (size_t ix = 0; ix < ply.size(); ++ix) {
    plyPoints_copy.push_back(make_pair(ix, ply[ix]));
  }
  cout << "finished making copy" << endl;

  size_t remaining_points = ply.size();

// Input : X = {x1, x2, ..., xn} , a set of 3D points (SimplePly ply)
//   P the number of planes to find (int nPlanes)
//   T the point−plane distance threshold (double threshold)
//   R the number of RANSAC trials (int nTrials)
// for p = 1 to P :
  for(int p = 1; p <= nPlanes; p++) {
//  bestPlane = {0 , 0}
    vector<PlyPoint> bestPlane;
    
//  bestPoints = {}
    vector<size_t> bestPoints;

    // size_t k = (sqrt(1 - w^nPlanes)) / (w^nPlanes);
    // std::cout << "k: " << k << std::endl;

    size_t inliers;
    double k;
//  for r = 1 to R :
    do {
      /*
//    S = {s1, s2, s3} = 3 points at random from X
      PlyPoint s1 = find_point(ply, null_colour);
      PlyPoint s2 = find_point(ply, null_colour);
      PlyPoint s3 = find_point(ply, null_colour);
      */
      PlyPoint s1 = plyPoints_copy[rand() % plyPoints_copy.size()].second;
      PlyPoint s2 = plyPoints_copy[rand() % plyPoints_copy.size()].second;
      PlyPoint s3 = plyPoints_copy[rand() % plyPoints_copy.size()].second;
      PlyPoint S[] = {s1, s2, s3};

//    thisPlane = {s1, crossProduct(s3-s1, s2-s1)}
      Eigen::Vector3d t1 = s3.location-s1.location;
      Eigen::Vector3d t2 = s2.location-s1.location;
      //Eigen::Vector3d n = t1.cross(t2);
      PlyPoint n;
      n.location = crossProduct(t1, t2);

      vector<PlyPoint> thisPlane;
      thisPlane.push_back(s1);
      thisPlane.push_back(n);

//    thisPoints = {}
      vector<size_t> thisPoints;
      
//    for xi in X :
      //for (size_t xi = 0; xi <= ply.size(); ++xi) {
      for (size_t xi = 0; xi <= plyPoints_copy.size(); ++xi) {
        // works in the same way as removing the points that have already been set to a plane.
        /** BAD DESIGN IDEA!!!! means that it's still iterating over all the points even when 
            they've been set to a plane. Making it take FOR EVER to run, realised to late meaning
            would have had to redo most of my RANSAC implentation so I ran with it.
        */
        //if(ply[xi].colour == null_colour) {

    //    if ( distance ( thisPlane , xi ) < T ) :
          //if (distance(thisPlane, plyP[xi]) < threshold) {

          if (distance(thisPlane, plyPoints_copy[xi].second) < threshold) {
    //      thisPoints = thisPoints + xi
            thisPoints.push_back(plyPoints_copy[xi].first);

            //cout << "plyPoints_copy before: " << plyPoints_copy.size() << endl;
            plyPoints_copy.push_back(plyPoints_copy[xi]);
            //cout << "plyPoints_copy during: " << plyPoints_copy.size() << endl;
            plyPoints_copy.erase(plyPoints_copy.begin() + xi);
            //cout << "plyPoints_copy after: " << plyPoints_copy.size() << endl;            
          }
          cout << "current Points: " << xi << endl;
        //}
      }
      cout << "thisPoints: " << thisPoints.size() << endl;
//    if | thisPoints | > | bestPoints | :
      if (thisPoints.size() > bestPoints.size()){
//      bestPlane = thisPlane
        bestPlane = thisPlane;
//      bestPoints = thisPoints
        bestPoints = thisPoints;
      }
      
      cout << "bestPoints: " << bestPoints.size() << endl;

      inliers = bestPoints.size();
      //std::cout << "inliers: " << inliers << std::endl;
      double w = (double)inliers / (double)remaining_points;
      k = sqrt(1 - pow(w, 3.0)) / pow(w, 3.0);
      //std::cout << "k: " << k << "\n" << std::endl;
      thisPoints.clear();
    } while (k > 45);

    for(unsigned point = 0; point < bestPoints.size(); point++) {
      size_t colourIx = p % colours.size(); // May need to recycle colours
      ply[bestPoints[point]].colour = colours[colourIx];
    }
    // for (size_t xi = 0; xi <= plyPoints_copy.size(); ++xi) {
    //   for(unsigned point = 0; point < bestPoints.size(); point++) {
    //     if (plyPoints_copy[xi].first == bestPoints[point]) {
    //       plyPoints_copy.erase(plyPoints_copy.begin() + xi);
    //       bestPoints.erase(bestPoints.begin() + point);
    //     }
    //   }
    // }

    plyPoints_copy = remove_points(plyPoints_copy, bestPoints.size());

    remaining_points = plyPoints_copy.size();
//  output bestPlane
    cout << "Plane"<< p << ": {" << endl;
    cout << "[" << bestPlane[0].location << "], " << endl;
    cout << "[" << bestPlane[1].location <<"]}"<< endl;

    cout << bestPoints.size() << endl;

    bestPoints.clear();
    bestPlane.clear();
//  X = X − bestPoints

  } // End RANSAC.

  // size_t planeSize = ply.size()/nPlanes;
  // for (size_t ix = 0; ix < ply.size(); ++ix) {
  //   size_t planeIx = ix / planeSize;
  //   size_t colourIx = planeIx % colours.size(); // May need to recycle colours
  //   ply[ix].colour = colours[colourIx];
  // }

  // Write the resulting (re-coloured) point cloud to a PLY file.
  cout << "Writing PLY data to " << argv[2] << endl;
  if (!ply.write(argv[2])) {
    cout << "Could not write PLY data to file " << argv[2] << endl;
    return -2;
  }
  return 0;
}

double distance(vector<PlyPoint> plane, PlyPoint p) {
  // d = unit(plane[1].location) . (plane[0].location - p.location)
  double result;
  Eigen::Vector3d unit_n = plane[1].location.normalized();
  result = unit_n.dot(plane[0].location - p.location);
  return fabs(result);
}

Eigen::Vector3d crossProduct(Eigen::Vector3d p1, Eigen::Vector3d p2) {
  Eigen::Vector3d result;
  result[0] = p1[1]*p2[2] - p1[2]*p2[1];
  result[1] = p1[2]*p2[0] - p1[0]*p2[2];
  result[2] = p1[0]*p2[1] - p1[1]*p2[0];
  return result;
}

/** Finds point that has still yet to have its colour set.*/
// PlyPoint find_point(SimplePly ply, Eigen::Vector3i null_colour) {
//   PlyPoint result;
//   // Not the best way to do this as, if all the points in the dataset have been set
//   // to a plane it could result in a never ending loop.

//   // (COME BACK TO).
//   while (true) {
//     PlyPoint temp = ply[rand() % ply.size()];
//     if (temp.colour == null_colour) return temp;
//   }
// }
vector<pair<size_t, PlyPoint> > remove_points(vector<pair<size_t, PlyPoint> > plyPoints, size_t points) {
  size_t count = points;
  while (count > 0) {
    plyPoints.pop_back();
    count--;
  }
  return plyPoints;
}
