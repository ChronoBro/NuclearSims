#ifndef _detector
#define _detector

#include "math.h"
#include "ring.h"
#include "ttt.h"
#include "random.h"

using namespace std;

struct location
{
  //double x,y,z;
  //double r,theta,phi;
  double steradian;

  double fx,fy,fz; //front and back stuff;
  double fr,fthetaMin,fphiMin, fthetaMax, fphiMax;
  double fdtheta, fdphi;
  
  double bx,by,bz;
  double br,bthetaMin,bphiMin, bthetaMax, bphiMax;
  double bdtheta, bdphi;
  
};

class detector
{
 public:

  detector(float, float,int); //distance to target along z-axis, thickness
  detector(float,float,float,float,int); //same as above but with x,y,z in lab frame specified (will add functionality later)
  ~detector();

  int event(float,float, int); //input theta and phi of particle, returns 1 if clean detection, 0 otherwise, last variable sets x,y,z or r,theta,phi check (0\1)
  int event(float,float,float,float,float,float, int); //same as above but gives velocity and t=0 position of decay
  
  int setGeometry(int,int,float,float,float,float,int); //# of x segments, # of y segments, xmin, xmax, ymin, ymax, detector type (0 or 1 for ring/rect)

  int setGeometry(int,int,float,float);
  int setGeometryFront(int,int,float,float);
  int setGeometryBack(int,int,float,float);
  int setGeometryFront(int,int,float,float,float,float,int);
  int setGeometryBack(int,int,float,float,float,float,int);
  
  int Nx;
  int Ny;
  float thetaHit;
  float phiHit;
  int segmentXhit;
  int segmentYhit;
  int didSetGeo;
  
 private:

  static CRandom ran;
  int coordinate = 0;
  //int is_hit;
  //float dist;
  float xD;
  float yD;
  float zD;
  float thickness;
  float theta;
  float phi;
  location ** Location;
  float dr;
  
};

#endif
