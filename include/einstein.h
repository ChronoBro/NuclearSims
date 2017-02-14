#ifndef einstein_
#define einstein_
#include <cmath>
#include <iostream>

class CEinstein
{
 public:
  float const c;
  float const nMass;
  CEinstein();
  void AddVelocities(float*, float*, float, float*);
  void FindCenterOfMass(float* , float, float*, float);
  float GetVelocity(float,float);
  void TwoBody(float,float,float,float&,float&);
  float theta;
  float phi;
  float energy;
  float velocity;
  float vcm[3];
  float velocitycm;
  float Erel;
};
#endif 
