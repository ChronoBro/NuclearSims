#include "twoFrag.h"


twoFrag::twoFrag(frag * part10, frag *part20, float Erel0, float gamma0)
{
  part1 = part10;
  part2 = part20;
  total_mass = part1->mass + part2->mass;
  reduced_mass = part1->mass*part2->mass/total_mass;
  Erel = Erel0;
  gamma = gamma0;
  Dvelocity = sqrt(2.*Erel/reduced_mass)*.9784;
  Vpart1 = Dvelocity*part2->mass/total_mass;
  Vpart2 = Dvelocity - Vpart1;
  Einstein = new einstein();
  Einstein->TwoBody(Erel,part1->mass,part2->mass,Vpart1,Vpart2);
}
//***********************************************************
void twoFrag::SampleWidth()
{
  float ErelSample = -1.;
  for (;;)
    { 
    ErelSample = part1->Array.Tower[0].Tel[0].gRandom.BreitWigner(Erel,gamma); 
    if (ErelSample > 0. && ErelSample < Erel + 10.*gamma) break;
    }

  Dvelocity = sqrt(2.*ErelSample/reduced_mass)*.9784;
  Vpart1 = Dvelocity*part2->mass/total_mass;
  Vpart2 = Dvelocity - Vpart1;
  Einstein->TwoBody(ErelSample,part1->mass,part2->mass,Vpart1,Vpart2);

}
//**********************************************************
float twoFrag::GetErel()
{
  /*
  float deltaVV = 0.;
  for (int i=0;i<3;i++) 
   {
     deltaV[i] = part1->v[i] - part2->v[i];
     deltaVV += pow(deltaV[i],2);
   }
   Recon.Dvelocity = sqrt(deltaVV);
   Recon.Erel = 0.5*reduced_mass*deltaVV/pow(.9784,2);
   return Recon.Erel;
  */
  Einstein->FindCenterOfMass(part1->v,part1->mass,part2->v,part2->mass);
  return Einstein->Erel;
}
//*******************************************************
float twoFrag::GetErelRecon()
{

  /*
  float deltaVV = 0.;
  float total_mom[3];
  float MM = 0.;
  for (int i=0;i<3;i++) 
   {
     deltaV[i] = part1->Recon.v[i] - part2->Recon.v[i];
     deltaVV += pow(deltaV[i],2);
     total_mom[i] = part1->Recon.v[i]*part1->mass + 
       part2->Recon.v[i]*part2->mass;
     MM += pow(total_mom[i],2);
   }
   Recon.Dvelocity = sqrt(deltaVV);
   Recon.Erel = 0.5*reduced_mass*deltaVV/pow(.9784,2);


   float M = sqrt(MM);
   Recon.plfTheta = acos(total_mom[2]/M);
   Recon.plfPhi = atan2(total_mom[1],total_mom[0]);
   Recon.plfVel = M/total_mass;
  */
  Einstein->FindCenterOfMass(part1->Recon.v,part1->mass,
  part2->Recon.v,part2->mass);
  Recon.plfTheta = Einstein->theta;
  Recon.plfPhi = Einstein->phi;
  Recon.plfVel = Einstein->velocity;

  Einstein->AddVelocities(part1->Recon.v,Einstein->vcm,part1->mass,
 part1->Recon.vcm);

  part1->Recon.velocitycm = 0.;
  for (int i=0;i<3;i++) part1->Recon.velocitycm += pow(part1->Recon.vcm[i],2);
  part1->Recon.velocitycm = sqrt(part1->Recon.velocitycm);

  float dot = 0.;
  for (int i=0;i<3;i++) dot += Einstein->vcm[i]*part1->Recon.vcm[i];
  Recon.thetaEmission = acos(-dot/Einstein->velocitycm/part1->Recon.velocitycm);
  Recon.plfVel = Einstein->velocitycm;

   return Einstein->Erel;

}
//**********************************************************
//checks to see that the two "detected" fragments are on top of each other
// i.e, in same Csi, or same si strips
int twoFrag::OnTopOf()
{
  int onTop = 0;
  if (part1->Array.hitTower != part2->Array.hitTower) return onTop;
  int itower = part1->Array.hitTower;
  if (part1->Array.Tower[itower].hitTele != part2->Array.Tower[itower].hitTele)
    return onTop;



  // at this stage both particles are in same telescope
  int itele = part1->Array.Tower[itower].hitTele;

  //these two telescopes have prolems preventing them from 
  //analysing two fragments in same telescope
  if (itower == 0 && itele == 3) return 1;
  if (itower == 1 && itele == 1) return 1;


  if (part1->Array.Tower[itower].Tel[itele].ICsI ==
      part2->Array.Tower[itower].Tel[itele].ICsI) return 1;

  if (part1->Array.Tower[itower].Tel[itele].ixStrip ==
      part2->Array.Tower[itower].Tel[itele].ixStrip) return 1;

  if (part1->Array.Tower[itower].Tel[itele].iyStrip ==
      part2->Array.Tower[itower].Tel[itele].iyStrip) return 1;

  return 0;


}
