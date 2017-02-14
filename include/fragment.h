#ifndef _fragment
#define _fragment

#include <cmath>

using namespace std;

/**
 *!\brief information about a single fragment 
 *
 */

class fragment
{
 public:
  float mass;  //!< fragment mass in amu
  float Z;     //!< fragment atomic number
  float velocity;  //!< fragment velocity in cm/ns
  float energy;    //!< fragment kinetic energy in MeV
  float theta;     //!<fragment theta angle in radians
  float phi;       //!<fragment phi angle in radians
  float v[3];      //!<lab velocity vector (cm/ns)
  float vcm[3];    //!<velocity vector in center-of-mass frame (cm/ns) 

  float pp;        //!< momentum
  float p[3];       //!< momentum vector
  float pcm[3];     //!< momentum vector in com

  fragment(){};
  fragment(float,float);
  void init(float,float);
  ~fragment();
  static float const pi;  //!< 3.14159.....
};

#endif
