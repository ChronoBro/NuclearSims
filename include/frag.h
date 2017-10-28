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
#include "detector.h"
#include "MathGeoLib.h"

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

  void setPosition(float *);
  void SetVelocity(float*,float*);
  void AddVelocity(float*);
  void AddVelocityRel(float*);
  float Eloss(float); // energy loss in target
  float Egain(float); // corrects for energy loss in target
  void MultiScat(float);
  bool targetInteraction(float,float);
  float p[3];
  void propagate(double tau); //put time in ns. Remember speed is in cm/ns but geometry is in mm!!!!!
  
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

  void setTrajectory(float x, float y, float z);
  void setRecon();
  void setRecon(float d, float measuredD, float start);
  Ray trajectory;
  double tofToRKE(double d, double TOF);
  double RKEToTOF(float RKE, float d);

  float ringDist;
  
};

#endif
