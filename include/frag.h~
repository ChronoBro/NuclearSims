#ifndef _frag
#define _frag

#include "array.h"
#include "loss.h"
#include <string>
#include "mScat.h"
#include "random.h"
#include "frame.h"
#include "fragment.h"
#include "ring.h"
#include "ring2.h"
#include "plane_det.h"

using namespace std;

/**
 *!\brief information about a single fragment and its interection with the detector
 *
 */
class CFrag : public fragment
{
 public:
  static CRandom ran;
  static float const pi; //!< 3.14159....

  float CsI_res;   

  CFrag(float,float,string,float,float);
  ~CFrag();
  CArray * Array;
  CArray * shadowArray;
  int hit();
  int hit3(float xtarget, float ytarget, float zbreakup, float straggle);
  int hit4(float xtarget, float ytarget, float zbreakup);
  int hit5(float xtarget, float ytarget, float zbreakup);
  int hit(float,float);
  int hitShadow(float,float);
  int is_hit;
  float ran_strag;
  float ran_phi;
  float x;
  float x_CsI;
  float y;
  float y_CsI;
  float rho;
  float rho_CsI;
  float theta_prime;
  float theta_prime_CsI;
  float phi_prime;
  float phi_prime_CsI;

  void SetVelocity(float*,float*);
  void AddVelocity(float*);
  void AddVelocityRel(float*);
  float Eloss(float); // energy loss in target
  float Egain(float); // corrects for energy loss in target
  void MultiScat(float);
  bool targetInteraction(float,float);

    CLoss *loss;
    CMultScat *multScat;

    CFrame *real;      //<!real particles energy, velocity, etc
    CFrame *recon;    //<!reconstructed properties

    bool alphaHole(int,int,int,int);
    bool protonHole(int,int,int,int);

    ring * Ring;
    ring * RingCsI;
    ring * Inactive;
    plane_det * Plane;

};

#endif
