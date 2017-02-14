#include "einstein.h"

using namespace std;

CEinstein::CEinstein():c(30.),nMass(931.478)
{

}
//*******************
//adds velocities in special relativity
// vel - velocity vector of object in some refence frame
// Vreference is ythe velocity of the refence frame relative to final frame
//Vfinal is velocity in final reference frame 

void CEinstein::AddVelocities(float* vel, float *Vreference, float mass, 
         float* Vfinal)

{
  // find parallel and perpendicular velocities to Vreference
  float Vpara[3];
  float Vperp[3];
  float dot = 0.;
  float VVreference = 0.;
  for (int i=0;i<3;i++) 
    {
      dot+= vel[i]*Vreference[i];
      VVreference += pow(Vreference[i],2);
    }

  for (int i=0;i<3;i++) 
    {
     Vpara[i] = dot/VVreference*Vreference[i];
     Vperp[i] = vel[i] - Vpara[i];
    }


  //now transform each
  float VfinalPara[3];
  float VfinalPerp[3];

  float cc = pow(c,2);
  float bb = sqrt(1.-VVreference/cc);

  float vv = 0.;
  for (int i=0;i<3;i++)
    {
      VfinalPara[i] = (Vpara[i] + Vreference[i])/(1.+dot/cc);
      VfinalPerp[i] = Vperp[i]*bb/(1.+dot/cc);
      Vfinal[i] = VfinalPara[i] + VfinalPerp[i];
      vv += pow(Vfinal[i],2);
    }
  velocity = sqrt(vv);
  theta = acos(Vfinal[2]/velocity);
  phi = atan2(Vfinal[1],Vfinal[0]);
  float pc = mass*nMass*velocity/c/sqrt(1.-vv/cc);
  float Etot = pc*c/velocity;
  energy = Etot - mass*nMass;
}

//************************************************************************
void CEinstein::FindCenterOfMass(float* v1, float mass1, float* v2 , 
    float mass2)
{
  /*
  cout << "v1R= " << v1[0] << " " << v1[1] << " " << v1[2] << endl;
  cout << "V2R="  << v2[0] << " " << v2[1] << " " << v2[2] << endl;
  */

  float vv1=0.;
  float vv2 = 0.;
  for (int i=0;i<3;i++)
    {
      vv1 += pow(v1[i],2);
      vv2 += pow(v2[i],2);
    }
  vv1 = sqrt(vv1);
  vv2 = sqrt(vv2);
  float p1 = mass1*nMass*vv1/c/sqrt(1.-pow(vv1/c,2));
  float p2 = mass2*nMass*vv2/c/sqrt(1.-pow(vv2/c,2));


  //define momentum vectors
  float pc1[3];
  float pc2[3];
  float pcTot[3];
  float pp = 0.;
  for (int i=0;i<3;i++)
    {
      pc1[i] = v1[i]*p1/vv1;
      pc2[i] = v2[i]*p2/vv2;
      pcTot[i] = pc1[i] + pc2[i];
      pp += pow(pcTot[i],2);
    }
  float Etot =  sqrt(pp + pow((mass1+mass2)*nMass,2));
  pp = sqrt(pp);
  float vvcm = pp/Etot*c;




  for (int i=0;i<3;i++) vcm[i] = -pcTot[i]/pp*vvcm;

  //cout << "cm= " << vcm[0] << " " << vcm[1] << " " << vcm[2] << endl;


  //transform back to center of frame
   float v1cm[3], v2cm[3];
   AddVelocities(v1,vcm,mass1,v1cm);
   float energy1 = energy;
   AddVelocities(v2,vcm,mass2,v2cm);
   float energy2 = energy;
   /*
   cout << "v1cm=" << v1cm[0] << " " << v1cm[1] << " " << v1cm[2] << endl;
   cout << "v2cm=" << v2cm[0] << " " << v2cm[1] << " " << v2cm[2] << endl;
   */
   Erel = energy1 + energy2;

   // define center of mass quantities
   velocitycm = vvcm;
   theta = acos(-vcm[2]/vvcm);
   phi = atan2(vcm[1],vcm[0]);
}
//**************************************************************
float CEinstein::GetVelocity(float energy, float mass)
{
  float Etot = energy + mass*nMass;
  float pc = sqrt(pow(Etot,2) - pow(mass*nMass,2));

  //cout << "energy= " << energy << endl;
  //cout << "vel = " << pc*c/Etot << endl;
  return pc*c/Etot;

}
//********************************************************
void CEinstein::TwoBody(float Erel, float mass1, float mass2, float&v1, 
float&v2)
{
  float Etot = (mass1+mass2)*nMass + Erel;
  double vv1 = v1;
  double EEtot = Etot;
  double cc = c;
  double nnMass = nMass;
  double mmass1 = mass1;
  double mmass2 = mass2;
  for (;;)
    {
      double pc = mmass1*nnMass*vv1/cc/sqrt(1.-pow(vv1/cc,2));
      double dpc_dv1 = pc/vv1 + 0.5*pc/(1.-pow(vv1/cc,2))*2.*vv1/cc/cc;
      double e1 = pc/vv1*cc;
      double de1_dv1 = dpc_dv1/vv1*cc - pc/vv1/vv1*cc;
      double e2 = sqrt(pow(pc,2)+pow(mmass2*nnMass,2));
      double de2_dv1 = pc*dpc_dv1/e2;
      v2 = pc/e2*c;
      double esum = e1 + e2;
      if (abs(esum-EEtot) < 1.e-6) break;
      double desum_dv1 = de1_dv1 + de2_dv1;
      vv1 -= (esum-Etot)/desum_dv1;
    }

}
