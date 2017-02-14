#ifndef scat_
#define scat_

#include <cmath>
#include <iostream>

using namespace std;

/**
 *!\brief Multiple scattering of a particle in the target
 */

class CMultScat
{
 public:
 
  float a; //!<screening parameter in fermi's
  float totalThick; //!<target thickness in mg/cm2
  float totalTau;
  float factor;
  static float const a0;
  static float const pi; //!< 3.14159......

  int Zpro;  //!< atomic number of projectile
  int Ztar;  //!< atomic number of target

  CMultScat(int,int,float); 
  float thetaRMS(float,float);
};

#endif
