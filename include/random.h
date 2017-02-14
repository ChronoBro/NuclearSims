#ifndef random_
#define random_
#include "TRandom3.h"
#include <cmath>

class CRandom
{
 protected:
  static TRandom3 random;
 public:

  static double Rndm();
  static float Gaus(float,float);
  static float expDecayTime(float);
  static float BreitWigner(float,float);
};


#endif
