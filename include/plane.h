#ifndef _plane
#define _plane

#include "math.h"
#include <iostream>


class plane
{
 public:

  plane();
  plane(double *n0, double * p0);
  ~plane();
  int setPlaneCoordinates(double * normX0, double * normY0);
  //void setBoundary(double * p10, double * p20, double * p30, double * p40);
  void findNormal(double *p10, double *p20, double * p30);
  double * labCoordinate(double xprime, double yprime); //give it plane coordinates
  void setPlaneCenter(double * p0);
  void setNormal(double * n0);
  void moveCenter(double, double, double);
  
 private:

  double n[3]; // define a plane by normal vector to plane, and point on the plane;
  double p[3];
  double normX[3];
  double normY[3];
  double labVector[3];
  int didSet;
  void checkSet();
  //double p1[3];
  //double p2[3];
  //double p3[3];
  //double p4[3];
  
  
};



#endif
