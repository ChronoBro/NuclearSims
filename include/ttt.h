#ifndef _ttt
#define _ttt
#include "random.h"
#include <iostream>
using namespace std;

class ttt
{
 public:
  static CRandom ran;
  ttt(float dist,float theta_D, float height, float width, float offset, int N_v, int N_h); //What are the units for these? How do I specify them?
  int hit (float theta, float phi);

  float thetaHit; //Reported theta value
  float phiHit; // Reported phi value

  float dist; //distance from target to center of detector in y-z plane, also target is assummed to be at 0 for cartesian map
  float theta_D; //Angle between distance vector and beam axis in y-z plane, z-axis is along beam, x-axis is along detector height and y-axis is defined accordingly for regular spherical coordinates 
  float height; //detector height
  float width; // detector width
  int N_v; //number of vertical strips
  int N_h; // number of horizontal strips
  int hitx; // horizontal section hit and x is vertical coordinate!
  int hity; // vertical section hit and y is horizontal coordinate!
  float a; //parameter1 for theta-phi map of detector plane
  float b; //parameter2 for theta-phi map of detector plane
  float c; //parameter3 for theta-phi map of detector plane (quadratic formula)
  float chi; 
  float x_min;
  float x_max;
  float x_1;//will be used for y extrema
  float x_2;
  float y_max;
  float y_min;
  float dy;
  float dx;
  float dl; //length transform for plane equation z=f(y)=-tan(theta_D)y + chi
  float rho;
  float offset; //offset of detector center from defined y-z plane
  float x;
  float y;
  float z;
  float rho1;
  float rho2;
  int CsIstrip1;
  int CsIstrip2;
};
#endif
