#ifndef frame_
#define frame_
#include <cmath>


/**
 *!\brief Relativistic and Non-Relativistic stuff
 *
 * give velocity, energies, transforms to new reference frame
 * using either non-relativistic or Relativistics equations
 */

class CFrame
{
 public:
  float energy;  //!< fragments energy in MeV
  float velocity; //!< fragments velocity in cm/ns 
  float pcTot;  //!< total momentum*c of fragment
  float v[3];  //!< velocity vector of fragment
  float pc[3]; //!< momentum*c vector of fragment MeV
  float mass;  //!< rest mass of fragment in MeV
  float A;
  float gamma; //!< relativistic gamma factor 
  float theta; //!< theta nagle of fragment in radians
  float phi;  //!< phi angle
  float totEnergy; //!< rest mass plas kinetic energy
  static float const m0; //!< nucleon mass in MeV
  static float const c; //!< speed of light (cm/ns)
  static float const konstant;
  static float const pi; //!< 3.14159....
 
  float getVelocity();
  float getVelocityNewton();
  float getVelocityRel();
  float getEnergy();
  float getEnergyNewton();
  float getEnergyRel();
  void getAngle();
  void transformVelocity(float*);
  void transformVelocityNewton(float*);
  void transformVelocityRel(float*);
  void getVelocityFromMom();
  void getVelocityFromMomNewton();
  void getVelocityFromMomRel();

  CFrame(float);
  static bool  einstein;
};


#endif
