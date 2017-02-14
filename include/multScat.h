#ifndef multScat_
#define multScat_
#include <cmath>

class multScat
{
 public:
  float const pi;
  float const alpha; //fine structure constant
  float const e2; //Coulomb const

  float Lrad;
  float LdashRad;
  float f;
  float totalThick; //atoms per cm2
  float Ztarget;
  float factor;

  float thetaRMS(float,float,float);
  multScat(float,float);

 };
#endif
