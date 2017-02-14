#include "frag.h"
#include "fragment.h"
#include <valarray>
#include "random.h"

struct prop
{
  float Erel;
  float Dvelocity;
  float plfTheta;
  float plfPhi;
  float plfVel;
  float thetaEmission;
};

class CDecay
{
 public:
  static CRandom ran;
  frag *part[5];
  frag *partCM[5];
  static float const pi;

  CDecay(frag*,frag*,frag*,frag*,frag*);
  ~CDecay();
  float GetErel();
  void ModeThru6Be();
  void ModeMicroCanonical();
  void ModeMicroCanonicalBe();
  void micro(int,fragment*,float,float);
  void getEk6Be();
  float getEk3body(frag*,frag*,frag*);
  void getEk6BeSec();
  float getEk3bodySec(frag*,frag*,frag*);
  int OnTopOf();
  prop Recon;

  float aaRel;
  float ppRel;
  float aapRelMin;
  float aapRelMax;
  float pEnergyMin;
  float pEnergyMax;
  float ppThetaRel;
  static float const EkTot8C;
  static float const gamma8C;
  float Ek6Be[6];
};
