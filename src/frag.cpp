#include "frag.h"



float const CFrag::pi=acos(-1.);
CRandom CFrag::ran;

/**
 * Constructor
 \param Z0 is the atomic number of fragment
 \param mass0 is mass of fragment in amu
 \param filename specifies the name of the file with energy loss information
 \param CsI_res0 fractional energy resolution of this fragment in CsI
 \param thickness is target thickness in mg/cm2 
*/

CFrag::CFrag(float Z0,float mass0, string filename, float CsI_res0,
	     float thickness):fragment(Z0,mass0)
{
  loss = new CLoss(filename,mass);
  CsI_res = CsI_res0;
  //CsI_res = 0;
  //be9 target  
  float thick = thickness/1000./9.*6.02e23; // atoms/ cem2
  multScat = new CMultScat((int)Z,4,thick);

  float fromCenter = 37.2;//35.;
  float radius = 63.; //50.;
  Array = new CArray(radius,fromCenter,6.42,7.6);
  shadowArray = new CArray(radius+2.,fromCenter,7.4,8.2);
  real = new CFrame(mass);
  recon = new CFrame(mass);

  float dist = 405.;
  //float dist = 370.;
  Ring = new ring(dist,7.5,62.5,128,128,0); //S4
  Inactive = new ring(dist,6.5,7.5,1,1,0);//S4 inactive area
  // RingCsI = new ring(365.,7.5,62.5,128,128); //S4 back of CsI
  RingCsI = new ring(410 + 39.,7.5,62.5,128,128,1); //p range;
  //RingCsI1 = new ring(332.,7.5,62.5,128,128);
  //RingCsI = new ring(315.,7.5,62.5,128,128); // 0 thickness

  Plane = new plane_det(445.,20.,20.,0.95);


}
//*********************************************************
/**
*Destructor
*/
CFrag::~CFrag()
{
  delete Array;
  delete real;
  delete recon;
}

//**********************************************************
  /**
   * logical function determines if a particle is detected
   */
int CFrag::hit()
{
  is_hit = Array->hit(real->theta,real->phi,(float)0.,(float)0.) ;

  if (is_hit) 
    {

      recon->theta = Array->thetaRecon;
      recon->phi = Array->phiRecon;

      recon->energy = real->energy + sqrt(real->energy)*CsI_res*
       ran.Gaus(0.,1.);

      //recon->energy = real->energy + real->energy*CsI_res*
      // ran.Gaus(0.,1.);




     recon->getVelocity();

    }

  return is_hit;
}
//**********************************************************
  /**
   * logical function determines if a particle is detected
   */
int CFrag::hit3(float xtarget,float ytarget, float zbreakup, float straggle)
{ 

  if (zbreakup > Ring->dist)
    {
      return 0;
    }
  rho = (Ring->dist - zbreakup)/cos(real->theta);
  rho_CsI = (RingCsI->dist - zbreakup)/cos(real->theta);

  ran_strag = ran.Gaus(0.,1.)*straggle;
  ran_phi = 2*3.14159*ran.Rndm();

  x = rho*sin(real->theta)*cos(real->phi) + xtarget;
  x_CsI = rho_CsI*sin(real->theta)*cos(real->phi) + xtarget + ran_strag*cos(ran_phi);

  y = rho*sin(real->theta)*sin(real->phi) + ytarget;
  y_CsI = rho_CsI*sin(real->theta)*sin(real->phi) + ytarget + ran_strag*sin(ran_phi);

  theta_prime = atan2(sqrt(pow(x,2.) + pow(y,2.)),Ring->dist);
  phi_prime = atan2(y,x);// + acos(-1);

  theta_prime_CsI = atan2(sqrt(pow(x_CsI,2.) + pow(y_CsI,2.)),410+39*cos(real->theta));
  phi_prime_CsI = atan2(y_CsI,x_CsI);// + acos(-1);
 //  if (phi_prime < 0. || phi_prime > 2.*acos(-1))
//     {
//       cout << "phi_prime = " << phi_prime << endl;
//       cout << "phi_prime_CsI = " << phi_prime_CsI << endl;
//       cout << endl;
//     }
//  cout << "phi here is " <<  phi_prime << endl;
  //is_hit = Ring->hit(theta_prime,phi_prime);
  is_hit = Ring->hit(real->theta, real->phi);
  //  cout << endl;
  //  cout << "now going to CsI " << endl;
  //  cout << "phi_CsI is " << phi_prime_CsI << endl;
  int dum1 = RingCsI->hit(theta_prime_CsI,phi_prime_CsI);
  //  cout << endl;


  //cout << "Csi res = " << CsI_res << endl;
  if (is_hit == 1 && dum1 == 1) 
    {

      //cout << "culprit? Ring->thetaHit = " << Ring->thetaHit << endl;
      recon->theta = Ring->thetaHit;
      recon->phi = Ring->phiHit;

      recon->energy = real->energy + sqrt(real->energy)*CsI_res*
	ran.Gaus(0.,1.);

      //cout << "culprit? recon->energy = " << recon->energy << endl;
      //recon->energy = real->energy + real->energy*CsI_res*
      // ran.Gaus(0.,1.);


     recon->getVelocity();
     return 1;
    }

  else 
    {
      return 0;
    }
}

// For heavy fragment in planar detector
  
int CFrag::hit4(float xtarget,float ytarget, float zbreakup)
{
  is_hit = Plane->hit(real->theta,real->phi,xtarget,ytarget,zbreakup) ;

  if (is_hit) 
    {

      recon->theta = Plane->thetaHit;
      recon->phi = Plane->phiHit;

      recon->energy = real->energy;

      // recon->energy = real->energy + real->energy*0.02*(2*ran.Rndm() - 1);

      //recon->energy = real->energy + sqrt(real->energy)*CsI_res*
      // ran.Gaus(0.,1.);

      //recon->energy = real->energy + real->energy*CsI_res*
      // ran.Gaus(0.,1.);




     recon->getVelocity();

    }

  return is_hit;
}


//For heavy fragment in inactive area
int CFrag::hit5(float xtarget, float ytarget, float zbreakup)
{
  if (zbreakup > Ring->dist)
    {
      return 0;
    }
  rho = (Inactive->dist - zbreakup)/cos(real->theta);
  x = rho*sin(real->theta)*cos(real->phi) + xtarget;
  y = rho*sin(real->theta)*sin(real->phi) + ytarget;
  theta_prime = atan2(sqrt(pow(x,2.) + pow(y,2.)),Ring->dist);
  phi_prime = atan2(y,x)+acos(-1);
  is_hit = Inactive->hit(theta_prime,phi_prime);
   return is_hit;
}

//*******************************************************
int CFrag::hitShadow(float xtarget , float ytarget)
{
  is_hit = shadowArray->hit(real->theta,real->phi,xtarget,ytarget) ;

  return is_hit;
}


//**********************************************************
  /**
   * logical function determines if a particle is detected
\param xTarget is in-plane location of interection in target form nomimal center
\param yTarget is out-of-plane location (cm)
  */
int CFrag::hit(float xTarget, float yTarget)
{
  is_hit = Array->hit(real->theta,real->phi,xTarget, yTarget) ;


  if (!is_hit) return is_hit;
  int itower = Array->hitTower;
  int itele = Array->Tower[itower]->hitTele;


  if (itower == 1 && itele == 0 && Array->Tower[itower]->Tele[itele]->ICsI == 2) return 0;


  if (itower == 2 && itele == 2 && Array->Tower[itower]->Tele[itele]->ICsI == 2) return 0;


  if (real->A == 1.)is_hit = protonHole(itower,itele,Array->Tower[itower]->Tele[itele]->ixStrip,31-Array->Tower[itower]->Tele[itele]->iyStrip);
  else if (real->A >3. )is_hit = alphaHole(itower,itele,Array->Tower[itower]->Tele[itele]->ixStrip,31-Array->Tower[itower]->Tele[itele]->iyStrip);


  if (is_hit)
    {
      
      recon->theta = Array->thetaRecon;
      recon->phi = Array->phiRecon;

      recon->energy = real->energy + sqrt(real->energy)*CsI_res*
       ran.Gaus(0.,1.);

      //recon->energy = real->energy + real->energy*CsI_res*
      // ran.Gaus(0.,1.);
     
     /*
      recon->theta = real->theta;
      recon->phi = real->phi;
      
      recon->energy = real->energy;
     */
     recon->getVelocity();

    }

  return is_hit;
}
//******************************************************************
  /**
   * Add a velocity vector to the fragments velocity vector.
   * Used to transform between reference frames
   */
void CFrag::AddVelocity(float *Vplf)
{
  real->transformVelocity(Vplf);
}
//****************************************************************
  /** 
   * returns the energy after the fragment has exited the target
\param thick is the thickness of target that the particle has to traverse (mg/cm2)
  */
float CFrag::Eloss(float thick)
{
  real->energy = loss->getEout(real->energy,thick);
  return real->energy;
}
//*******************************************************************
  /**
   * corrects energy of a detected particle for the energy loss
   * in the target.
\param thick is the thickness of target material though which the particle passed (mg/cm2)
  */
float CFrag::Egain(float thick)
{
  if (thick > 0.)
    recon->energy = loss->getEin(recon->energy,thick/cos(recon->theta));

  recon->getVelocity();

   return recon->energy;
}
//***********************************************
//include multiple scattering
  /**
   * Monte Carlo choice of small angle scattering due to passage through the target
\param fractionalThick is the fractional thick of the target through which the particle passed
  */
void CFrag::MultiScat(float fractionalThick)
{

  float thetaRMS = multScat->thetaRMS(real->energy,fractionalThick);
  float sigma = thetaRMS/sqrt(2.);
  //cout << "thetaRMS= " << thetaRMS << endl;
  float deltaTheta = sqrt(2.)*sigma*sqrt(-log(ran.Rndm()));
  //cout << "deltaTheta= " << deltaTheta << endl;
  float deltaPhi = 2.*pi*ran.Rndm();
  //cout << "delta Phi= " << deltaPhi << endl;

  float x = sin(deltaTheta)*cos(deltaPhi);
  float y = sin(deltaTheta)*sin(deltaPhi);
  float z = cos(deltaTheta);



  //rotate in z-x plane by theta
  float xx = x*cos(real->theta) + z*sin(real->theta);
  float yy = y;
  float zz = z*cos(real->theta) - x*sin(real->theta);


  //rotate in x-y plane
  float xxx = xx*cos(real->phi) - yy*sin(real->phi);
  float yyy = yy*cos(real->phi) + xx*sin(real->phi);
  float zzz = zz;


  float thetaNew = acos(zzz);
  float phiNew = atan2(yyy,xxx);


  real->theta = thetaNew;
  real->phi = phiNew;
}
//*********************
  /**
   * accounts for multiscattering and energy loss in the target
   \param dthick is thickness of target though the particle must pass (mg/cm2)
\param thickness is total target thickness (mg/cm2)
   */
bool CFrag::targetInteraction(float dthick, float thickness)
{
  bool stopped = 0;
  if (dthick == 0.) return stopped;
  float thick = dthick/cos(real->theta);
  Eloss(thick);
  if (real->energy <= 0.) 
    {
      stopped = 1;
      return stopped;
    }
  MultiScat(thick/thickness);
  return stopped;
}
//**************************************************
bool CFrag::protonHole(int itower, int itele, int ifront, int iback)
{


  if (itower == 0 && itele == 1) //det0
    {
      if (iback == 20) return false;
      if (iback == 30) return false;
      if (ifront == 20) return false;
      if (ifront == 30)  return false;

      float x = 1.;

      if (iback == 0) x *= .87;
      if (iback == 15) x *= .74;
      if (iback == 16) x *= .77;
      if (iback == 24) x *= .68;
      if (iback == 31) x *= .88;
      if (ifront == 0) x *= .87;
      if (ifront == 15) x *= .77;
      if (ifront == 16) x *= .77;
      if (ifront == 24) x *= .68;
      if (ifront == 31) x *= .86;
      if (ran.Rndm() > x) return false;
    }
  else if (itower == 0 && itele == 0)//det1
    {
      float x = 1.;
      if (iback == 0) x *= .84;
      if (iback == 15) x *= .69;
      if (iback == 16) x *= .88;
      if (iback == 18) x *= .95;
      if (iback == 31) x *= .86;
      if (ifront == 0) x *= .87;
      if (ifront == 14) x *= .85;
      if (ifront == 15) x *= .68;
      if (ifront == 31) x *= .89;
      if (ran.Rndm() > x) return false;
    }

  else if (itower == 1 && itele == 2)//det2
    {
      float x = 1.;
      if (iback == 0) x *= .98;
      if (iback == 15) x *= .80;
      if (iback == 16) x *= .77;
      if (iback == 17) x*= .87;
      if (iback == 31) x*= .92;
      if (ifront == 0) x*= .91;
      if (ifront == 15) x*= .72;
      if (ifront == 16) x*= .88;
      if (ifront == 31) x*=.94;
      if (ran.Rndm() > x) return false;
    }
  else if (itower == 1 && itele == 1) //det3
    {
      float x = 1.;
      if (iback == 0) x *= .88;
      if (iback == 15) x*= .67;
      if (iback == 16) x*= .87;
      if (iback == 31) x*= .86;
      if (ifront == 0) x*= .88;
      if (ifront == 15) x*= .68;
      if (ifront == 16) x*= .88;
      if (ifront == 31) x*= .86;
      if (ran.Rndm() > x) return false;
    }
  else if (itower == 1 && itele == 0) //det4
    {
      if (iback == 19) return false;
      if (ifront == 19)  return false;
      float x = 1.;
      if (iback == 0) x *= .90;
      if (iback == 15) x*= .66;
      if (iback == 31) x*= .91;
      if (ifront == 0) x*=.91;
      if (ifront == 15) x*=.64;
      if (ifront == 31) x*= .89;
      if (ran.Rndm() > x) return false;

    }
  else if (itower == 2 && itele == 3) //det5
    {
      float x = 1.;
      if (iback == 0) x*= .87;
      if (iback == 16) x*=.71;
      if (iback == 17) x*=.86;
      if (iback == 31) x*=.81;
      if (ifront == 0) x*=.91;
      if (ifront == 14) x*=.96;
      if (ifront == 15) x*=.76;
      if (ifront == 16) x*=.69;
      if (ifront == 31) x*=.91;
      if (ran.Rndm() > x) return false;
    }
  else if (itower == 2 && itele == 2) //det 6
    {
      float x = 1.;
      if (iback == 0) x*= .86;
      if (iback == 16) x*=.70;
      if (iback == 31) x*=.72;
      if (ifront == 0) x*=.83;
      if (ifront == 15) x*=.81;
      if (ifront == 16) x*=.71;
      if (ifront == 31) x*=.85;
      if (ran.Rndm() > x) return false;
    }
  else if (itower == 2 && itele == 1) //det 7
    {
      if (iback == 4) return false;
      float x = 1.;
      if (iback == 0) x*= .86;
      if (iback == 15) x*=.66;
      if (iback == 16) x*=.85;
      if (iback == 31) x*=.66;
      if (ifront == 0) x*=.86;
      if (ifront == 15) x*=.71;
      if (ifront == 16) x*=.81;
      if (ifront == 31) x*=.83;
      if (ran.Rndm() > x) return false;
    }
  else if (itower == 2 && itele == 0) //det 8
    {
      float x = 1.;
      if (iback == 0) x *= .83;
      if (iback == 14) x*=.89;
      if (iback == 15) x*=.66;
      if (iback == 31) x*=.87;
      if (ifront == 0) x*=.86;
      if (iback == 15) x*=.71;
      if (iback == 16) x*=.81;
      if (iback == 31) x*=.83;
      if (ran.Rndm() > x) return false;
    }
  else if (itower == 3 && itele == 2) //det 9
    {
      if (iback == 27) return false;
      float x = 1.;
      if (iback == 0) x *= .93;
      if (iback == 15) x*= .78;
      if (iback == 31) x*=.85;
      if (ifront == 0)x*=.92;
      if (ifront == 31)x*=.95;
      if (ran.Rndm() > x) return false;
    }
  else if (itower == 3 && itele == 1) //det 10
    {
      float x = 1.;
      if (iback == 0) x *= .89;
      if (iback == 15) x*= .68;
      if (iback == 16) x*= .81;
      if (iback == 31) x*= .86;
      if (ifront == 0) x*=.89;
      if (ifront == 15) x*=.90;
      if (ifront == 16) x*=.62;
      if (ifront == 17) x*=.94;
      if (ifront == 31) x*=.86;
      if (ran.Rndm() > x) return false;
    }
  else if (itower == 3 && itele == 0) //det 11
    {
      float x = 1.;
      if (iback == 0) x*= .87;
      if (iback == 14) x*=.93;
      if (iback == 15) x*=.63;
      if (iback == 16) x*=.94;
      if (iback == 31) x*=.89;
      if (ifront == 0) x*=.92;
      if (ifront == 15) x*=.67;
      if (ifront == 16) x*= .90;
      if (ifront == 31) x*=.91;
      if (ran.Rndm() > x) return false;
    }
  else if (itower == 4 && itele == 1) //det 11
    {
      float x = 1.;
      if (iback == 0) x*= .89;
      if (iback == 15) x*= .88;
      if (iback == 16) x*=.75;
      if (iback == 31) x*=.90;
      if (ifront == 0) x*=.92;
      if (ifront == 16) x*=.71;
      if (ifront == 17) x*=.83;
      if (ifront == 31) x*=.84;
      if (ran.Rndm() > x) return false;
    }
  else if (itower == 4 && itele == 0) //det 11
    {
      if (iback == 31) return false;
      float x = 1.;
      if (iback == 0) x*=.86;
      if (iback == 9) x*= .42;
      if (iback == 15) x*=.67;
      if (iback == 16) x*=.90;
      if (ifront == 0) x*=.93;
      if (ifront == 16) x*=.68;
      if (ifront == 17) x*=.83;
      if (ifront == 31) x*=.87;
      if (ran.Rndm() > x) return false;
    }

  return true;
}
//********************************************
bool CFrag::alphaHole(int itower, int itele, int ifront, int iback)
{



  if (itower == 0 && itele == 1) //det0
    {
      if (iback == 20) return false;
      if (iback == 30) return false;
      float x = 1.;

      if (iback == 0)  x *= .86;
      if (iback == 15) x *= .79;
      if (iback == 16) x *= .91;
      if (iback == 31) x *= .82;
      if (ifront == 0) x *= .75;
      if (ifront == 15) x*= .68;
      if (ifront == 31) x*= .84;
      if (ran.Rndm() > x) return false;

    }
  else if (itower == 0 && itele == 0)//det1
    {
      float x = 1.;
      if (iback == 0) x *= .87;
      if (iback == 15) x *=.81;
      if (iback == 31) x *= .86;
      if (ifront == 0) x *= .75;
      if (ifront == 15) x*= .69;
      if (ifront == 31) x *= .84;
      if (ran.Rndm() > x) return false;
    }
  else if (itower == 1 && itele == 2)//det2
    {
      float x = 1.;
      if (iback == 16) x *= .76;
      if (ifront == 15) x *= .71;
      if (ifront == 31) x *= .92;
      if (ran.Rndm() > x) return false;
    }
  else if (itower == 1 && itele == 1) //det3
    {
      float x = 1.;
      if (iback == 0) x *= .86;
      if (iback == 15) x *= .74;
      if (iback == 16) x *= .95;
      if (iback == 31) x *= .85;
      if (ifront == 0) x *= .86;
      if (ifront == 15) x *= .74;
      if (ifront == 16) x *= .96;
      if (ifront == 31) x *= .84;
      if (ran.Rndm() > x) return false;
    }
  else if (itower == 1 && itele == 0) //det4
    {
      float x = 1.;
      if (iback == 19) return false;
      if (ifront == 19)  return false;
      if (iback == 0) x *= .86;
      if (iback == 15) x *= .72;
      if (iback == 31) x *= .90;
      if (ifront == 0) x *= .86;
      if (ifront == 15) x *= .72;
      if (ifront == 31) x *= .88;
      if (ran.Rndm() > x) return false;
    }
  else if (itower == 2 && itele == 3) //det5
    {
      float x = 1.;
      if (iback == 16) x *= .62;
      if (ifront == 0) x *= .9;
      if (ifront == 15) x *= .82;
      if (ifront == 16) x *= .72;
      if (ifront == 31) x *= .87;
      if (ran.Rndm() > x) return false;
    }
  else if (itower == 2 && itele == 2) //det6
    {
      float x  = 1.;
      if (iback == 0) x *= .88;
      if (iback == 16) x *= .75;
      if (iback == 31) x *= .85;
      if (ifront == 0) x *= .81;
      if (ifront == 15) x *= .9;
      if (ifront == 16) x *= .81;
      if (ifront == 31) x *= .82;
      if (ran.Rndm() > x) return false;
    }
  else if (itower == 2 && itele == 1) //det7
    {
      if (iback == 14) return false;
      float x = 1.;
      if (iback == 0) x *= .80;
      if (iback == 15) x *= .76;
      if (iback == 31) x *= .87;
      if (ifront == 0) x *= .81;
      if (ifront == 15) x *= .76;
      if (ifront == 16) x *= .9;
      if (ifront == 31) x *= .8;
      if (ran.Rndm() > x) return false;
    }
  else if (itower == 2 && itele == 0) //det8
    {
      float x = 1.;
      if (iback == 0) x *= .82;
      if (iback == 15) x *= .78;
      if (iback == 31) x *= .91;
      if (ifront == 0) x *= .86;
      if (ifront == 15) x *= .76;
      if (ifront == 16) x *= .90;
      if (ifront == 31) x *= .83;
      if (ran.Rndm() > x) return false;
    }
  else if (itower == 3 && itele == 2) //det9
    {
      if (iback == 27) return false;
      float x = 1.;
      if (iback == 0) x *= .91;
      if (iback == 15) x *= .84;
      if (iback == 31) x *= .86;
      if (ifront == 0) x *= .87;
      if (ifront == 0) x *= .90;
      if (ran.Rndm() > x) return false;
    }
  else if (itower == 3 && itele == 1) //det10
    {
      float x = 1.;
      if (iback == 0) x *= .85;
      if (iback == 15) x *= .74;
      if (iback == 16) x *= .92;
      if (iback == 31) x *= .84;
      if (ifront == 0) x *= .86;
      if (ifront == 16) x *= .71;
      if (ifront == 31) x *= .79;
      if (ran.Rndm() > x) return false;

    }
  else if (itower == 3 && itele == 0) //det11
    {
      float x = 1.;
      if (iback == 0) x*= .83;
      if (iback == 15) x *= .73;
      if (iback == 31) x *= .87;
      if (ifront == 0) x *= .86;
      if (ifront == 15) x *= .78;
      if (ifront == 31) x *= .91;
      if (ran.Rndm() > x) return false;
    }
  else if (itower == 4 && itele == 1) //det12
    {
      float x = 1.;
      if (iback == 0) x *= .88;
      if (iback == 31) x *= .90;
      if (ifront == 0) x *= .89;
      if (ifront == 16) x *= .73;
      if (ifront == 31) x *= .84;
      if (ran.Rndm() > x) return false;
    }
  else if (itower == 4 && itele == 0) //det13
    {
      if (iback == 31)  return false;
      float x = 1.;
      if (iback == 0) x *= .82;
      if (iback == 15) x *= .74;
      if (ifront == 0) x *= .94;
      if (ifront == 16) x *= .74;
      if (ifront == 31) x *= .82;
      if (ran.Rndm() > x) return false;
    }





  return true;
}
