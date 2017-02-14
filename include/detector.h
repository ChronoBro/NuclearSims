#ifndef _detector
#define _detector

#include "math.h"
#include "ring.h"
#include "loss.h"
#include "ttt.h"

using namespace std;

/* struct location */
/* { */
/*   double x,y,z; */
/*   double r,theta,phi; */
/*   double deltatheta; */
/*   double deltaphi; */
/*   double steradian; */
/* }; */

/* struct teleP */
/* { */
/*   location Location[64][48]; */
/* }; */

class detector
{
 public:

  detector(float, float); //distance to target along z-axis, thickness
  detector(float,float,float,float); //same as above but with x,y,z in lab frame specified (will add functionality later)
  ~detector();
    
  int hit(float,float); //input theta and phi of particle, returns 1 if clean detection, 0 otherwise
  
  int setGeometry(int,int,float,float,float,float,int); //# of x segments, # of y segments, xmin, xmax, ymin, ymax, detector type (0 or 1 for ring/rect)

  int setGeometry(int,int,float,float);
 private:

  int is_hit;
  float dist;
  float x;
  float y;
  float z;
  float thickness;
  float theta;
  float phi;

  float ** xy;
  
};

#endif
