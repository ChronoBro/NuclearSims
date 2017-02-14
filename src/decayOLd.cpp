#include "decay.h"

CRandom CDecay::ran;
float const CDecay::pi = acos(-1.);
float const CDecay::EkTot8C = 3.51;
float const CDecay::gamma8C = .23;
CDecay::CDecay(frag * part10, frag *part20, frag* part30, frag* part40, 
	     frag* part50)
{
  part[0] = part10;
  part[1] = part20;
  part[2] = part30;
  part[3] = part40;
  part[4] = part50;

  string blank("");
  for (int i=0;i<5;i++) 
    {
     partCM[i] = new frag(part[i]->Z,part[i]->mass,blank,0.,0.);
     partCM[i]->mass = part[i]->mass;
    }
}
//*********************************************************
CDecay::~CDecay()
{
  for (int i=0;i<4;i++) delete partCM[i];
}
//**********************************************************
float CDecay::GetErel()
{



  valarray <float> vcm(3);
  for (int i=0;i<3;i++) 
   {
     for (int j=0;j<5;j++)
       vcm[i] += part[j]->Recon.v[i]*part[j]->mass;
   }
  vcm /= 8.;  
  float vv = 0.;
  for (int i=0;i<3;i++) vv += pow(vcm[i],2);
  Recon.plfVel = sqrt(vv);
  Recon.plfTheta = acos(vcm[2]/Recon.plfVel);
  Recon.plfPhi = atan2(vcm[1],vcm[0]);

  Recon.Erel = 0.;
   for (int j=0;j<5;j++)
     {
       partCM[j]->velocity = 0.;
       for (int i=0;i<3;i++)
         {
          partCM[j]->v[i] = part[j]->Recon.v[i] - vcm[i];
          partCM[j]->velocity += pow(partCM[j]->v[i],2);
         }
       partCM[j]->energy = part[j]->mass/2.*partCM[j]->velocity/pow(.9784,2);
       partCM[j]->velocity = sqrt(partCM[j]->velocity);
       Recon.Erel += partCM[j]->energy;
     }

   return Recon.Erel;

}


//*************************************************************
void CDecay::ModeMicroCanonical()
{


  valarray <float> vcm(3);
  for (int i=0;i<5;i++)
    {
      part[i]->velocity = ran.Gaus(0.,1.)/sqrt(part[i]->mass); //velocity
      part[i]->pp = part[i]->velocity*part[i]->mass; //momentum
      part[i]->theta = acos(1.-2.*ran.Rndm());
      part[i]->phi = 2.*pi*ran.Rndm();

      part[i]->v[0] = part[i]->velocity*sin(part[i]->theta)*cos(part[i]->phi);
      part[i]->v[1] = part[i]->velocity*sin(part[i]->theta)*sin(part[i]->phi);
      part[i]->v[2] = part[i]->velocity*cos(part[i]->theta);
      
      for (int j=0;j<3;j++) vcm[j] += part[i]->v[j]*part[i]->mass;
    }

  vcm /= 8.;

  float testTotal= 0.;
  for (int i=0;i<5;i++)
    {
      part[i]->velocity = 0.;
      for (int j=0;j<3;j++)
	{
	  part[i]->v[j] -= vcm[j];
          part[i]->velocity += pow(part[i]->v[j],2);
         }
      part[i]->energy = part[i]->mass/2.*part[i]->velocity/pow(.9784,2);
      part[i]->velocity = sqrt(part[i]->velocity);
      testTotal += part[i]->energy;
    }
  float ratio = sqrt(EkTot8C/testTotal);
  for (int i=0;i<5;i++)
    {
      part[i]->velocity *= ratio;
      for (int j=0;j<3;j++) part[i]->v[j] *= ratio;
    }


}

//************************************************************
int CDecay::OnTopOf()
{

  for (int i=0;i<4;i++)
    for (int j=i+1;j<5;j++)
      {
	if (part[i]->Array->hitTower == part[j]->Array->hitTower)
	  {
	    int itower = part[i]->Array->hitTower;
            if (part[i]->Array->Tower[itower]->hitTele 
                 == part[j]->Array->Tower[itower]->hitTele) 
	      {
		int itele = part[i]->Array->Tower[itower]->hitTele;
                int icsi1 = part[i]->Array->Tower[itower]->Tele[itele]->ICsI;
                int icsi2 = part[j]->Array->Tower[itower]->Tele[itele]->ICsI;
		if (icsi1 == icsi2) return 1;
                int ix1 = part[i]->Array->Tower[itower]->Tele[itele]->ixStrip;
                int ix2 = part[j]->Array->Tower[itower]->Tele[itele]->ixStrip;
		if (ix1 == ix2) return 1;
                int iy1 = part[i]->Array->Tower[itower]->Tele[itele]->iyStrip;
                int iy2 = part[j]->Array->Tower[itower]->Tele[itele]->iyStrip;
		if (iy1 == iy2) return 1;
	      }

	  }
      }

  return 0;
}
  //*****************************
    void CDecay::micro(int N,fragment *Frag,float Ektot,float massTot)
{
  valarray <float> vcm(3);
  for (int i=0;i<N;i++)
    {

      Frag[i].v[0] = ran.Gaus(0.,1.)/Frag[i].mass;
      Frag[i].v[1] = ran.Gaus(0.,1.)/Frag[i].mass;
      Frag[i].v[2] = ran.Gaus(0.,1.)/Frag[i].mass;
      
      for (int j=0;j<3;j++) vcm[j] += Frag[i].v[j]*Frag[i].mass;
    }

  vcm /= massTot;

  float testTotal= 0.;
  for (int i=0;i<N;i++)
    {
      Frag[i].velocity = 0.;
      for (int j=0;j<3;j++)
	{
	  Frag[i].v[j] -= vcm[j];
          Frag[i].velocity += pow(Frag[i].v[j],2);
         }
      Frag[i].energy = Frag[i].mass/2.*Frag[i].velocity/pow(.9784,2);
      Frag[i].velocity = sqrt(Frag[i].velocity);
      testTotal += Frag[i].energy;
    }
  float ratio = sqrt(Ektot/testTotal);
  for (int i=0;i<N;i++)
    {
      Frag[i].velocity *= ratio;
      for (int j=0;j<3;j++) Frag[i].v[j] *= ratio;
    }

}

//**********************************************************
//sequential decay throough first excited state of 9B
void CDecay::ModeThru6Be()
{

  //avergae total kinetic energy for 6Be decay
  float const EkTot6BeMean = 1.3711;
  // width 
  float const EkTot6BeGamma = .092;

  //choose total kinetic 
  float EkTot;
  for (;;)
    {
      EkTot = ran.BreitWigner(EkTot8C,gamma8C);
      if( EkTot > EkTot6BeMean + 2.*EkTot6BeGamma) break;
    }




  //first choose the total kinetic energy for 6Be decay
  float EkTot6Be;
  for(;;)
    {
      EkTot6Be = ran.BreitWigner(EkTot6BeMean,EkTot6BeGamma);
      if (EkTot6Be > 0.) break;
    }
    

  //first decay 6Be+2p
  float EkTotFirst = EkTot - EkTot6Be;

  fragment Frag[3];
  Frag[0].init(1.,1.);
  Frag[1].init(1.,1.);
  Frag[2].init(4.,6.);
  fragment recoil;
  recoil.init(4,6);

  micro(3,Frag,EkTotFirst,8.);
  for (int i=0;i<3;i++)
    {
      part[0]->v[i] = Frag[0].v[i];
      part[1]->v[i] = Frag[1].v[i];
      recoil.v[i] = Frag[2].v[i];
    } 

  //second decay
  Frag[2].init(4.,6.);
  micro(3,Frag,EkTot6Be,6.);
  for (int i=0;i<3;i++)
    {
      part[2]->v[i] = Frag[0].v[i] + recoil.v[i];
      part[3]->v[i] = Frag[1].v[i] + recoil.v[i];
      part[4]->v[i] = Frag[2].v[i] + recoil.v[i];
    } 

}
//*********************************************************
float CDecay::getEk3body(frag* frag1, frag* frag2, frag* frag3)
{
  float totMass = frag1->mass + frag2->mass + frag3->mass;

  float vcm[3] = {0.};
  float ek = 0;
  for (int i=0;i<3;i++)
    {
      vcm[i] += frag1->v[i]*frag1->mass;
      vcm[i] += frag2->v[i]*frag2->mass;
      vcm[i] += frag3->v[i]*frag3->mass;
      vcm[i] /= totMass;

      frag1->vcm[i] = frag1->v[i] - vcm[i];
      frag2->vcm[i] = frag2->v[i] - vcm[i];
      frag3->vcm[i] = frag3->v[i] - vcm[i];
      
      ek += 0.5*frag1->mass*pow(frag1->vcm[i]/.9784,2);
      ek += 0.5*frag2->mass*pow(frag2->vcm[i]/.9784,2);
      ek += 0.5*frag3->mass*pow(frag3->vcm[i]/.9784,2);
    }
  return ek;
}
//************************************************************
void CDecay::getEk6Be()
{

  Ek6Be[0] = getEk3body(part[0],part[1],part[4]);
  Ek6Be[1] = getEk3body(part[0],part[2],part[4]);
  Ek6Be[2] = getEk3body(part[0],part[3],part[4]);
  Ek6Be[3] = getEk3body(part[1],part[2],part[4]); 
  Ek6Be[4] = getEk3body(part[1],part[3],part[4]);
  Ek6Be[5] = getEk3body(part[2],part[3],part[4]);

}
//*********************************************************
float CDecay::getEk3bodySec(frag* frag1, frag* frag2, frag* frag3)
{
  float totMass = frag1->mass + frag2->mass + frag3->mass;

  float vcm[3] = {0.};
  float ek = 0;
  for (int i=0;i<3;i++)
    {
      vcm[i] += frag1->Recon.v[i]*frag1->mass;
      vcm[i] += frag2->Recon.v[i]*frag2->mass;
      vcm[i] += frag3->Recon.v[i]*frag3->mass;
      vcm[i] /= totMass;

      frag1->vcm[i] = frag1->Recon.v[i] - vcm[i];
      frag2->vcm[i] = frag2->Recon.v[i] - vcm[i];
      frag3->vcm[i] = frag3->Recon.v[i] - vcm[i];
      
      ek += 0.5*frag1->mass*pow(frag1->vcm[i]/.9784,2);
      ek += 0.5*frag2->mass*pow(frag2->vcm[i]/.9784,2);
      ek += 0.5*frag3->mass*pow(frag3->vcm[i]/.9784,2);
    }
  return ek;
}
//************************************************************
void CDecay::getEk6BeSec()
{

  Ek6Be[0] = getEk3bodySec(part[0],part[1],part[4]);
  Ek6Be[1] = getEk3bodySec(part[0],part[2],part[4]);
  Ek6Be[2] = getEk3bodySec(part[0],part[3],part[4]);
  Ek6Be[3] = getEk3bodySec(part[1],part[2],part[4]); 
  Ek6Be[4] = getEk3bodySec(part[1],part[3],part[4]);
  Ek6Be[5] = getEk3bodySec(part[2],part[3],part[4]);

}
