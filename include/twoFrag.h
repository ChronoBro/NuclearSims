#include "frag.h"


struct prop
{
  float Erel;
  float Dvelocity;
  float plfTheta;
  float plfPhi;
  float plfVel;
  float thetaEmission;
};

class twoFrag
{
 public:
  frag *part1;
  frag *part2;

  twoFrag(frag*,frag*,float,float);
  float GetErel();
  float GetErelRecon();
  void SampleWidth();
  int OnTopOf();



  float deltaV[3];
  float Dvelocity;
  float Erel;
  float gamma;
  float reduced_mass;
  float total_mass;
  float Vpart1;
  float Vpart2;

  prop Recon;
  einstein *Einstein;

};
