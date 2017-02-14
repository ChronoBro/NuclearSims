#include "frag.h"
#include "fragment.h"
#include <valarray>
#include "random.h"
#include "moscow.h"

struct prop
{
  float Erel;
  float Dvelocity;
  float plfTheta;
  float plfPhi;
  float plfVel;
  float thetaEmission;
};


/**
 *!\brief selects the veloity vectors of the secondary fragments
 */

class CDecay
{
 public:
  CFrame* frag6Be;
  float CosTheta_T;
  float x_T;
  float x_Y[2];
  float CosTheta_Y[2];

  void getJacobi(CFrame**,bool);
  void getJacobiPrimary();
  void getJacobiSecondary();
  


  bool einstein;
  static CRandom ran;
  CFrame *real[3]; 
  CFrame *recon[3];
  CFrame *plfRecon;
  CFrame *partCM[3];
  CFrag *frag[3];  //!< information about the decay fragments
  static float const pi; //!< 3.14159....

  CDecay(CFrag*,CFrag*,CFrag*,bool einstein0);
  ~CDecay();
  float getErelReal();
  float getErelRecon();
  float getErel(CFrame**);
  float getErelNewton(CFrame**);
  float getErelRel(CFrame**);

  void ModeMoscow();
  void ModeMicroCanonical();
  void ModeMicroCanonicalBe();
  void micro(int,CFrame**,float,float);
  void getEk6Be(CFrame**);
  float getEk3body(CFrame*,CFrame*,CFrame*);
  float getEk3bodyNewton(CFrame*,CFrame*,CFrame*);
  float getEk3bodyRel(CFrame*,CFrame*,CFrame*);
  float getEk3();
  void getEk6BeRecon();
  void getEk6BeReal();
  bool OnTopOf();
  bool OnTopOf3();
  bool leaveCsI();
  prop Recon;

  float  sumA;
  float ErelRecon; //!<reconstructed relative kinetic energy


  float aaRel;
  float ppRel[6];
  float aapRelMin;
  float aapRelMax;
  float pEnergyMin;
  float pEnergyMax;
  float ppThetaRel;
  static float const EkTot8B;
  static float const gamma8B;
  float Ek6Be[6];
  float aRatio;

  int Nsolution; //!< number of p-p pairs with correct 6Be energy
  int Isolution; //!< solution #

  moscow * Moscow;
};
