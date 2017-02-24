#include "plf.h"
using namespace std;

CRandom CPlf::ran;
float const CPlf::pi = acos(-1.);

/**
 * Constructor
 \param thickness is target thickness in mg/cm2
*/

CPlf::CPlf(float thickness)
{
  temp1 = .5*pi/180.; 
  temp2 = 1.4*pi/180.;
  frame = new CFrame(8.);
  float thick = thickness/1000./9.*6.02e23; // atoms/ cem2
  multScat =  new CMultScat(20,4);

  p[0] = 0.;
  p[1] = 0.;
  p[2] = 0.;
  
}
//********************************************
  /**
   * Destructor
   */
CPlf::~CPlf()
{
  delete frame;
  delete multScat;
}
//********************************************
  /**
   * Monte Carlo selection of the angle and velocity of the parent fragment
   */
void CPlf::GetPlf()
{
  // float vv;
  // if (ran.Rndm()>0.2) vv = ran.Gaus(11.01,.2);
  // else 
  // {
  //    float y1 = ran.Gaus(0.,1);
  //    float y2 = ran.Gaus(0.,1);
  //    float y = sqrt(pow(y1,2)+pow(y2,2));
  //    vv = 10.9 - y/3.;
  //  }

  // float phi = 2.*pi*ran.Rndm();
  // float theta;
  // if(ran.Rndm()>0.1)theta = sqrt(temp1)*sqrt(-log(ran.Rndm()+1.e-32));
  // else theta = sqrt(temp2)*sqrt(-log(ran.Rndm()+1.e-32));

  

  //  float v[3];
  //  v[0] = vv*sin(theta)*cos(phi);
  //  v[1] = vv*sin(theta)*sin(phi);
  //  v[2] = vv*cos(theta);

  //  for(int i=0;i<3;i++)frame->v[i] = v[i];
  //  frame->velocity= vv;
  //  frame->theta = theta;
  //  frame->phi=phi;

  float vv;
  float ptr;
  float pz;
  float pbeam;
  float pmag;

  //CMomDist momDist;
  
  ptr = 0.; //momDist.getTransMom(); could make these gaussian or something
  pz = 0.; //momDist.getLongMom();
  pbeam = 12508.1;
  
  float theta = 0.;
  theta = atan2(ptr, pz + pbeam);

  float phi = 0.;
  phi = 2.*pi*ran.Rndm();

  pmag = sqrt(pow(pz + pbeam,2.)+pow(ptr,2.));
  vv = pmag/sqrt(pow(pmag,2.)+pow(31964.9,2.));
  vv = 30.*vv;
  
  float v[3];
  v[0] = vv*sin(theta)*cos(phi);
  v[1] = vv*sin(theta)*sin(phi);
  v[2] = vv*cos(theta);
  
  for (int i=0;i<3;i++)frame->v[i] = v[i];
  frame->velocity = vv;
  frame->theta = theta;
  frame->phi = phi;

}

void CPlf::GetPlf(float EPA0, int mass, int dontChange)
{
  // float vv;
  // if (ran.Rndm()>0.2) vv = ran.Gaus(11.01,.2);
  // else 
  // {
  //    float y1 = ran.Gaus(0.,1);
  //    float y2 = ran.Gaus(0.,1);
  //    float y = sqrt(pow(y1,2)+pow(y2,2));
  //    vv = 10.9 - y/3.;
  //  }

  // float phi = 2.*pi*ran.Rndm();
  // float theta;
  // if(ran.Rndm()>0.1)theta = sqrt(temp1)*sqrt(-log(ran.Rndm()+1.e-32));
  // else theta = sqrt(temp2)*sqrt(-log(ran.Rndm()+1.e-32));

  

  //  float v[3];
  //  v[0] = vv*sin(theta)*cos(phi);
  //  v[1] = vv*sin(theta)*sin(phi);
  //  v[2] = vv*cos(theta);

  //  for(int i=0;i<3;i++)frame->v[i] = v[i];
  //  frame->velocity= vv;
  //  frame->theta = theta;
  //  frame->phi=phi;

  //below is isotropic breakup
  
  float vv;
  float ptr;
  float pz;
  float pbeam;
  float pmag;

  //CMomDist momDist;

  if(!dontChange)
    {
      ptr = momDist.getTransMom(); //could make these gaussian or something
      pz = momDist.getLongMom();
    }
  else
    {
      ptr=0;
      pz=0;
    }
  
  pbeam = sqrt(pow((EPA0+931.478)*mass,2)-pow(931.478*(float)mass,2));
  
  float theta = 0.;
  theta = atan2(ptr, pz + pbeam);

  float phi = 0.;
  phi = 2.*pi*ran.Rndm();

  pmag = sqrt(pow(pz + pbeam,2.)+pow(ptr,2.));
  vv = pmag/sqrt(pow(pmag,2.)+pow(931.494*(float)mass,2.));
  vv = 30.*vv;
  
  float v[3];
 
  if(!dontChange)
    {
      v[0] = vv*sin(theta)*cos(phi);
      v[1] = vv*sin(theta)*sin(phi);
      v[2] = vv*cos(theta);
      for (int i=0;i<3;i++)frame->v[i] = v[i];
    }
  else
    {
       for (int i=0;i<3;i++)frame->v[i] = frame->v[i]*vv/frame->velocity;
    }
  
  frame->velocity = vv;

  if(!dontChange)
    {
      frame->theta = theta;
      frame->phi = phi;
    }
  //cout << "pbeam = " << pbeam << endl;
  
}



//*************************************************
  /**
   * selects breakup angle of parent fragments (isotropic distribution)
   */
void CPlf::isotropic()
{
  theta_breakup = acos(1.-2.*ran.Rndm());
  phi_breakup = 2.*pi*ran.Rndm();

  unit_breakup[0] = sin(theta_breakup)*cos(phi_breakup);
  unit_breakup[1] = sin(theta_breakup)*sin(phi_breakup);
  unit_breakup[2] = cos(theta_breakup);

}
//*************************************************
  /**
   * Selects break angle of parent fragment.
   * the distribution is assumed to be isotropic in-plane i.e. 1/sin(theta)
   */
void CPlf::inplane()
{
  theta_breakup =pi*ran.Rndm();
  phi_breakup = 2.*pi*ran.Rndm();

  unit_breakup[0] = sin(theta_breakup)*cos(phi_breakup);
  unit_breakup[1] = sin(theta_breakup)*sin(phi_breakup);
  unit_breakup[2] = cos(theta_breakup);


  //rotate along axis of plf
  //first the x-z plane
  float xd = unit_breakup[0]*cos(theta) + unit_breakup[2]*sin(theta);
  float yd = unit_breakup[1];
  float zd = unit_breakup[2]*cos(theta) - unit_breakup[0]*sin(theta);

  //then in x-y plane
  float xdd = xd*cos(phi) + yd*sin(phi);
  float ydd = yd*cos(phi) - xd*sin(phi);
  float zdd = zd;

  unit_breakup[0] = xdd;
  unit_breakup[1] = ydd;
  unit_breakup[2] = zdd;


}
//***********************************************
//include multiple scattering
  /**
   * multiple scattering of parent fragment
   \param fractionalThick fractional thickness of the target through which
    the parent passes before it decays (mg.cm2)
   */
void CPlf::MultiScat(float fractionalThick)
{
  if (fractionalThick == 0.) return;
  frame->getEnergy();
  float thetaRMS = multScat->thetaRMS(frame->energy,fractionalThick);
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
  float xx = x*cos(frame->theta) + z*sin(frame->theta);
  float yy = y;
  float zz = z*cos(frame->theta) - x*sin(frame->theta);


  //rotate in x-y plane
  float xxx = xx*cos(frame->phi) - yy*sin(frame->phi);
  float yyy = yy*cos(frame->phi) + xx*sin(frame->phi);
  float zzz = zz;


  float thetaNew = acos(zzz);
  float phiNew = atan2(yyy,xxx);
  if (phiNew < 0.) phiNew += 2.*pi;

  frame->theta = thetaNew;
  frame->phi = phiNew;
  frame->v[0] = frame->velocity*cos(phiNew)*sin(thetaNew);
  frame->v[1] = frame->velocity*sin(phiNew)*sin(thetaNew);
  frame->v[2] = frame->velocity*cos(thetaNew);
  
}

void CPlf::propagate(double tau)
{
  p[0] = tau*v[0]*10. + p[0]; //velocities are in cm/ns
  p[1] = tau*v[1]*10. + p[1];
  p[2] = tau*v[2]*10. + p[2];

}

void CPlf::setBeamSpot(double targetSize)
{

  do
    {
      p[0] = (1.-2.*ran.Rndm())/2.*targetSize;
      p[1] = (1.-2.*ran.Rndm())/2.*targetSize;
    } while( sqrt(pow(p[0],2)+pow(p[1],2)) > targetSize/2.);
  
  
  
}

void CPlf::setLifetime(double decayConstant)
{

  TRandom eRan(0);
  float lifetime = (float)eRan.Exp(frame->gamma*decayConstant);
  propagate((double)lifetime);
  
}
