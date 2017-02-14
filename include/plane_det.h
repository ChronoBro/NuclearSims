#ifndef _plane_det
#define _plane_det
#include "random.h"
#include <iostream>
using namespace std;

class plane_det
{
 public:
  static CRandom ran;
  plane_det(float dist, float half_width, float half_height, float det_res);
  int hit(float theta, float phi, float xtarget, float ytarget, float zbreakup);

  float thetaHit;
  float phiHit;

  float z;
  float xmax;
  float ymax;
  float res;

  float rho;
  float x;
  float y;
  float det_x;
  float det_y;
};
#endif
