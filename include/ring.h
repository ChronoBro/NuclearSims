#ifndef _ring
#define _ring
#include "random.h"
#include <iostream>
using namespace std;

class ring
{
 public:
  static CRandom ran;
  ring(float dist,float rmin, float rmax, int Nring, int Npi, int place);
  int hit(float theta, float phi);

  int hitRing; // ring which was hit
  int hitPie;  // pie sector which was hit

  float thetaHit; //zenith angle in radians assigned to hit
  float phiHit; // azimuth angle in radians assigned to hit

  float rmin;  //minium radius of ring counter
  float rmax;  //maximum radius of ring counter
  float dist;  //distance of ring counter from target
  int Nring;   // number of annullar ring sectors
  int Npie;    //number of pie sectors
  int place;
  float theta_min;
  float theta_max;
  std::vector<float> theta_min_vec;
  std::vector<float> theta_max_vec;
  float dtheta;
  float dphi;
  float CsIstrip;
  int segmentation_ring;
  int CsIring;
  float dr;
  int useful;
  int please;
};
#endif
