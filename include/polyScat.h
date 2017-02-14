#include "multScat.h"

class polyScat
{
 public:
  multScat* Hydrogen;
  multScat* Carbon;
  float Zhydrogen;
  float Zcarbon;
  float Nhydrogen;
  float Ncarbon;

  polyScat(float);
  float thetaRMS(float,float,float);
};
