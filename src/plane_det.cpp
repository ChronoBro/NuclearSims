#include "plane_det.h"

CRandom plane_det::ran;

plane_det::plane_det(float dist, float half_width, float half_height, float det_res)
{
  z = dist;
  xmax = half_width;
  ymax = half_height;
  res = det_res;

}

int plane_det::hit(float theta, float phi, float xtarget, float ytarget, float zbreakup)
{
  rho = (z - zbreakup)/cos(theta);
  x = rho*sin(theta)*cos(phi)+xtarget;
  y = rho*sin(theta)*sin(phi)+ytarget;
 

  if (abs(x) > xmax)
    {
      return 0;
    }
  else if (abs(y) > ymax)
    {
      return 0;
    }
  else
    {
      det_x = x + res*ran.Gaus(0.,1.);
      det_y = y + res*ran.Gaus(0.,1.);
      thetaHit = atan2(sqrt(pow(det_x,2.) + pow(det_y,2.)),z);
      phiHit = atan2(det_y,det_x);
      return 1;
    }
}
