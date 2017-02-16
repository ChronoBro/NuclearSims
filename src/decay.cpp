#include "decay.h"

CRandom CDecay::ran;
float const CDecay::pi = acos(-1.);
float const CDecay::EkTot8B = 2.842; // energy released in 8B IAS decay to 6LI IAS
float const CDecay::gamma8B = 0.00075;

/**
 * Constructor
\param part10 is a pointer to particle 1
\param part20 is a pointer to particle 2
\param part30 is a pointer to particle 3

 */
CDecay::CDecay(CFrag * part10, CFrag *part20, CFrag* part30, 
	       bool einstein0)
{

  Nfrags = 3;
  frag6Be = new CFrame (6);

  einstein = einstein0;
  CFrame::einstein = einstein;

  frag[0] = part10;
  frag[1] = part20;
  frag[2] = part30;



  real[0] = part10->real;
  real[1] = part20->real;
  real[2] = part30->real;



  recon[0] = part10->recon;
  recon[1] = part20->recon;
  recon[2] = part30->recon;


  sumA = 0.;
  for (int i=0;i<3;i++) sumA += real[i]->A;

  plfRecon = new CFrame(sumA);


  for (int i=0;i<3;i++) 
    {
     partCM[i] = new CFrame(real[i]->A);
    }
 

  string name("6be_p2_new.dat");
  Moscow = new moscow(name);
}

CDecay::CDecay(CFrag * part10, CFrag *part20,bool einstein0)
{

  Nfrags = 2;
  //frag6Be = new CFrame (6);

  einstein = einstein0;
  CFrame::einstein = einstein;

  frag[0] = part10;
  frag[1] = part20;


  real[0] = part10->real;
  real[1] = part20->real;


  recon[0] = part10->recon;
  recon[1] = part20->recon;

  sumA = 0.;
  for (int i=0;i<Nfrags;i++) sumA += real[i]->A;

  plfRecon = new CFrame(sumA);


  for (int i=0;i<Nfrags;i++) 
    {
     partCM[i] = new CFrame(real[i]->A);
    }
 
}


CDecay::CDecay(CFrag * part10, CFrag *part20,bool einstein0, float vCM0)
{

  Nfrags = 2;
  vCM = vCM0;
  einstein = einstein0;
  CFrame::einstein = einstein;

  frag[0] = part10;
  frag[1] = part20;


  real[0] = part10->real;
  real[1] = part20->real;


  recon[0] = part10->recon;
  recon[1] = part20->recon;

  sumA = 0.;
  for (int i=0;i<2;i++) sumA += real[i]->A;

  plfRecon = new CFrame(sumA);


  for (int i=0;i<2;i++) 
    {
     partCM[i] = new CFrame(real[i]->A);
    }
 
  for (int i=0;i<180;i++)
    {
      float theta = ((float)i+.5)*3.14159/180.;
      //float dsigma = pow(sin(theta),5);  // m= 2
      //float dsigma =sin(theta)*pow((3.*pow(cos(theta),2)-1.),2); // m = 0
      float dsigma = (1.-.92*pow(cos(theta),2))*sin(theta);
      angle[i] = dsigma;
      if (i > 0) angle[i] += angle[i-1];
    }

  for (int i=0;i<180;i++) 
    {
     angle[i] /= angle[179];
     //cout << angle[i] << endl;
    }


  for (int i=0;i<360;i++)
    {
      float phi = ((float)i+.5)*3.14159/180.;
      float dsigma = 3052.06;
      dsigma -= 1.11417e5*(exp(-pow((phi - pi/2.)/(23.5031*pi/180.),2)/2.)+exp(-pow((phi - 3.*pi/2.)/(23.5031*pi/180.),2)/2.));
      dsigma += 1.11464e5*(exp(-pow((phi - pi/2.)/(23.2806*pi/180.),2)/2.)+exp(-pow((phi - 3.*pi/2.)/(23.2806*pi/180.),2)/2.));
      phiAngle[i] = dsigma;
      if (i > 0) phiAngle[i] += phiAngle[i-1];
    }

  for (int i=0;i<360;i++) 
    {
     phiAngle[i] /= phiAngle[359];

    }

}




//*********************************************************
  /**
   * Destructor
   */
CDecay::~CDecay()
{
  for (int i=0;i<3;i++) delete partCM[i];
  delete plfRecon;
  delete Moscow;
}
//**********************************************************
  /**
   * returns the reconstructed kinetical energy of the fragmens in their
   * center-of-mass frame using the real fragment velocities.
   */
float CDecay::getErelReal()
{
  return getErel(real);
}
//********************************************************
//*********************************************************
  /**
   * returns the reconstructed kinetical energy of the fragmens in their
   * center-of-mass frame using the reconstructed or detected
   *  fragment velocities. 
   */
float CDecay::getErelRecon()
{
  return getErel(recon);
}
//*********************************************************
float CDecay::getErel(CFrame** part)
{
  if (einstein) return getErelRel(part);
  else return getErelNewton(part);
}
//**********************************************************
  /** 
   * find the relative kinetic energy of the fragments in their
   * center-of-mass frame. Non-relativistic version
   \param part is a pointer to the fragments velocity vectors (real or reconstructed)
   */

float CDecay::getErelNewton(CFrame** part)
{

  for (int i=0;i<3;i++) 
   {
     plfRecon->v[i] = 0.;
     for (int j=0;j<3;j++)
           plfRecon->v[i] += part[j]->v[i]*part[j]->A;
     plfRecon->v[i] /= sumA;  
   }

  plfRecon->getEnergy();

  ErelRecon = 0.;

   for (int j=0;j<Nfrags;j++)
     {
       partCM[j]->velocity = 0.;
       for (int i=0;i<3;i++)
         {
          partCM[j]->v[i] = part[j]->v[i] - plfRecon->v[i];
          partCM[j]->velocity += pow(partCM[j]->v[i],2);
         }
       partCM[j]->energy = real[j]->A/2.*partCM[j]->velocity/pow(.9784,2);
       partCM[j]->velocity = sqrt(partCM[j]->velocity);
       ErelRecon += partCM[j]->energy;
     }

   return ErelRecon;

}

//**********************************************************
  /** 
   * find the relative kinetic energy of the fragments in their
   * center-of-mass frame. Relativistic version
   \param part is a pointer to the fragments velocity vectors (real or reconstructed)
   */
float CDecay::getErelRel(CFrame **part)
{


  for (int i=0;i<3;i++) 
   {
     plfRecon->pc[i] = 0.;
     for (int j=0;j<Nfrags;j++)plfRecon->pc[i] += part[j]->pc[i];
   }
  plfRecon->totEnergy = 0;
  for (int j=0;j<Nfrags;j++) plfRecon->totEnergy += part[j]->totEnergy;
  plfRecon->getVelocityFromMom();


  plfRecon->velocity = sqrt(pow(plfRecon->v[0],2)+pow(plfRecon->v[1],2)
		    +pow(plfRecon->v[2],2));
  plfRecon->theta = acos(plfRecon->v[2]/plfRecon->velocity);


  for (int j=0;j<Nfrags;j++)
    {
    for (int i=0;i<3;i++)
      {
	partCM[j]->v[i] = part[j]->v[i];
      }
    }
  //we want to "subtract" the velocity to get to the CoM frame
  float vReference[3];
  for (int i=0;i<3;i++)
    {
      vReference[i] = -1*plfRecon->v[i];
    }

  
  
  ErelRecon = 0.;
   for (int j=0;j<Nfrags;j++)
     {
       partCM[j]->transformVelocity(vReference);
       ErelRecon += partCM[j]->getEnergy();
     }


   return ErelRecon;


   //p-p relative
   int ii = 0;
   for (int i=0;i<1;i++)
     {
       for (int j=i+1;j<2;j++)
	 {
           for (int k=0;k<3;k++) 
             {
               plfRecon->pc[k] = part[i]->pc[k]+part[j]->pc[k];
	     }
          plfRecon->totEnergy = part[i]->totEnergy
                                + part[j]->totEnergy;
          plfRecon->getVelocityFromMom();
          for (int k=0;k<3;k++) 
	    {
             partCM[i]->v[k] = part[i]->v[k];   
             partCM[j]->v[k] = part[j]->v[k];   
	    }
             ppRel[ii] = 0.;

             partCM[i]->transformVelocity(plfRecon->v);
             ppRel[ii] += partCM[i]->getEnergy();

             partCM[j]->transformVelocity(plfRecon->v);
             ppRel[ii] += partCM[j]->getEnergy();
	     ii++;
	 }
     }


  for (int i=0;i<3;i++) 
   {
     plfRecon->pc[i] = 0.;
     for (int j=0;j<4;j++)plfRecon->pc[i] += part[j]->pc[i];
   }
  plfRecon->totEnergy = 0;
  for (int j=0;j<2;j++) plfRecon->totEnergy += part[j]->totEnergy;
  plfRecon->getVelocityFromMom();


  for (int j=0;j<2;j++)
    for (int i=0;i<3;i++) partCM[j]->v[i] = part[j]->v[i];

  float ErelReconA = 0.;
   for (int j=0;j<2;j++)
     {
       partCM[j]->transformVelocity(plfRecon->v);
       ErelReconA += partCM[j]->getEnergy();
     }

   aRatio = ErelReconA/ErelRecon;


}
//*************************************************************
  /**
   * The momentum of the five fragments are chosen randomly from the 
   * allowable phase space.
   */
void CDecay::ModeMicroCanonical()
{

  float EkTot;
//   for(;;)
//     {
//       EkTot = ran.BreitWigner(EkTot8B,gamma8B);
//          if( fabs(EkTot- EkTot8B) <  2.*gamma8B && EkTot > 0.) break;
//     }
  //EkTot = 2.2; // ground state
  //EkTot = 1.476;
  EkTot = 1.96; 
  
  micro(real,EkTot,sumA);
}

//************************************************************
  /**
   * checks to see if each of the five fragments go into separate Si strips
   * (both front and back) and separate CsI crystals 
   */
bool CDecay::OnTopOf()
{


  
  //for inner two detectors consider 3x3 array of CsI
  for (int i=0;i<Nfrags;i++)
    {
      if (frag[i]->Array->hitTower != 2) continue;
      if (frag[i]->Array->Tower[2]->hitTele == 1)
        frag[i]->Array->Tower[2]->Tele[1]->ICsI = 
	  frag[i]->Array->Tower[2]->Tele[1]->ICsI3; 
      if (frag[i]->Array->Tower[2]->hitTele == 2)
        frag[i]->Array->Tower[2]->Tele[2]->ICsI = 
	  frag[i]->Array->Tower[2]->Tele[2]->ICsI3; 

      if (frag[i]->shadowArray->hitTower != 2) continue;
      if (frag[i]->shadowArray->Tower[2]->hitTele == 1)
        frag[i]->shadowArray->Tower[2]->Tele[1]->ICsI = 
	  frag[i]->shadowArray->Tower[2]->Tele[1]->ICsI3; 
      if (frag[i]->shadowArray->Tower[2]->hitTele == 2)
        frag[i]->shadowArray->Tower[2]->Tele[2]->ICsI = 
	  frag[i]->shadowArray->Tower[2]->Tele[2]->ICsI3; 
    }
  

  for (int i=0;i<2;i++)
    for (int j=i+1;j<Nfrags;j++)
      {
	if (frag[i]->Array->hitTower == frag[j]->Array->hitTower)
	  {
	    int itower = frag[i]->Array->hitTower;
            if (frag[i]->Array->Tower[itower]->hitTele 
                 == frag[j]->Array->Tower[itower]->hitTele) 
	      {
		int itele = frag[i]->Array->Tower[itower]->hitTele;
                int icsi1 = frag[i]->Array->Tower[itower]->Tele[itele]->ICsI;
                int icsi2 = frag[j]->Array->Tower[itower]->Tele[itele]->ICsI;
		if (icsi1 == icsi2) return 1;
                int ix1 = frag[i]->Array->Tower[itower]->Tele[itele]->ixStrip;
                int ix2 = frag[j]->Array->Tower[itower]->Tele[itele]->ixStrip;
		if (ix1 == ix2) return 1;
                int iy1 = frag[i]->Array->Tower[itower]->Tele[itele]->iyStrip;
                int iy2 = frag[j]->Array->Tower[itower]->Tele[itele]->iyStrip;
		if (iy1 == iy2) return 1;
	      }

	  }
      }

  return 0;
}

bool CDecay::OnTopOf3()
{

  for(int i=0;i<2;i++)
    {
      float csistrip = frag[i]->Ring->CsIstrip;
      float csiring = frag[i]->Ring->CsIring;
      float csistrip1 = frag[i]->RingCsI->CsIstrip;
      float csiring1 = frag[i]->RingCsI->CsIring;
      if(csistrip!=csistrip1 || csiring!=csiring1) return true;
    }
  
  if (frag[0]->Ring->CsIstrip == frag[1]->Ring->CsIstrip && frag[0]->Ring->CsIring == frag[1]->Ring->CsIring) return true;
  // if (frag[1]->Ring->CsIstrip == frag[2]->Ring->CsIstrip && frag[1]->Ring->CsIring == frag[2]->Ring->CsIring) return true;
  //if (frag[0]->Ring->CsIstrip == frag[2]->Ring->CsIstrip && frag[0]->Ring->CsIring == frag[2]->Ring->CsIring) return true;

  if (frag[0]->Ring->hitRing == frag[1]->Ring->hitRing) return true;
  //if (frag[0]->Ring->hitRing == frag[2]->Ring->hitRing) return true;
  //if (frag[1]->Ring->hitRing == frag[2]->Ring->hitRing) return true;
  return false;


}




//**********************************************************************
  /**
   * checks to see if each of the four protons leaves its CsI
   */
bool CDecay::leaveCsI()
{
  for (int i=0;i<2;i++) //only protons
    {
      if (frag[i]->Array->hitTower != frag[i]->shadowArray->hitTower)
	{
	  cout << "shadow hit tower does not match" << endl;
          abort();
	}
      int itower = frag[i]->Array->hitTower;
      if (frag[i]->Array->Tower[itower]->hitTele 
                 != frag[i]->shadowArray->Tower[itower]->hitTele) 
	{
	  cout << "shadow hit detector does not match" << endl;
          abort();

	}
      int itele = frag[i]->Array->Tower[itower]->hitTele;
      if (frag[i]->Array->Tower[itower]->Tele[itele]->ICsI !=  
	  frag[i]->shadowArray->Tower[itower]->Tele[itele]->ICsI) return 1;     
    }
  return 0;
}
  //*****************************
    void CDecay::micro(CFrame **Frag,float Ektot,float massTot)
{
  valarray <float> vcm(3);
  for (int i=0;i<Nfrags;i++)
    {

      Frag[i]->v[0] = ran.Gaus(0.,1.)/Frag[i]->A;
      Frag[i]->v[1] = ran.Gaus(0.,1.)/Frag[i]->A;
      Frag[i]->v[2] = ran.Gaus(0.,1.)/Frag[i]->A;
      
      for (int j=0;j<Nfrags;j++) vcm[j] += Frag[i]->v[j]*Frag[i]->A;
    }

  vcm /= massTot;

  float testTotal= 0.;
  for (int i=0;i<Nfrags;i++)
    {
      Frag[i]->velocity = 0.;
      for (int j=0;j<3;j++)
	{
	  Frag[i]->v[j] -= vcm[j];
          Frag[i]->velocity += pow(Frag[i]->v[j],2);
         }
      Frag[i]->energy = Frag[i]->A/2.*Frag[i]->velocity/pow(.9784,2);
      Frag[i]->velocity = sqrt(Frag[i]->velocity);
      testTotal += Frag[i]->energy;
    }
  float ratio = sqrt(Ektot/testTotal);
  for (int i=0;i<Nfrags;i++)
    {
      Frag[i]->velocity *= ratio;
      for (int j=0;j<Nfrags;j++) Frag[i]->v[j] *= ratio;
    }

}

void CDecay::Mode(double ET, double Ex, double mass1, double mass2)
{

  //double ET_mean = 0.3111  ;//Q = -7.268
  //double Excitation = .2439;
  //double Excitation = 0.;
  //double Et;
  //for (;;)
  //  {
  //    Et = ran.BreitWigner(ET_mean,.1);
  //    if (Et > .1  && Et < 10.) break;
      
  //  }

  double Et = ET + Ex;

  double mu = mass1*mass2/(mass1+mass2);
  double vrel = sqrt(2.*Et/mu)*.9784;
  double v1 = vrel*mass2/(mass1+mass2);
  double v2 = vrel - v1;


  double theta = acos(1.-2.*ran.Rndm());
  double phi = 2.*acos(-1.)*ran.Rndm();
  real[0]->v[0] = v1*sin(theta)*cos(phi);
  real[0]->v[1] = v1*sin(theta)*sin(phi);
  real[0]->v[2] = v1*cos(theta);


  for (int i=0;i<3;i++) real[1]->v[i] = -v2/v1*real[0]->v[i];



}




//**********************************************************
  /**
   * Simulates the decay of 8C as two sequential prompt two-proton decays
   * passing through the 6Be ground state. The prompt decays are each 
   * treated as sampling available phase space.
   */
void CDecay::ModeMoscow()
{

  //avergae total kinetic energy for 6Be decay
  float const EkTot6Be = 1.3711;
  // width 
  float const gamma6Be = .092;

  //choose total kinetic 
  float EkTot;
  if (gamma6Be == 0.)EkTot= EkTot6Be;
  else
    {
     for (;;)
       {
         EkTot = ran.BreitWigner(EkTot6Be,gamma6Be);
         if( fabs(EkTot- EkTot6Be) <  2.*gamma6Be) break;
        }
    }


  CFrame ** Frag;
  Frag = new CFrame* [3];
  Frag[0] = new CFrame(1.);
  Frag[1] = new CFrame(1.);
  Frag[2] = new CFrame(4.);

  
  Moscow->getEvent(EkTot,Frag[2],Frag[1],Frag[0]);

  //look at correlations in Jacobi coordinates


  for (int i=0;i<3;i++)
    {
      real[0]->v[i] = Frag[0]->v[i];
      real[1]->v[i] = Frag[1]->v[i];
      real[2]->v[i] = Frag[2]->v[i];
    } 
  getJacobiPrimary();
}



//********************************************************************
float CDecay::getEk3body(CFrame* frag1, CFrame* frag2, CFrame* frag3)
{
  if (einstein) return getEk3bodyRel(frag1,frag2,frag3);
  else return getEk3bodyNewton(frag1,frag2,frag3);
}

//*********************************************************************
  /**
   * returns the relative kinetic energy of three fragments in their
   * center-of-mass frame. Non-relativistic version
   \param frag1 is a pointer to fragment 1 
   \param frag2 is a pointer to fragment 2 
   \param frag3 is a pointer to fragment 3 
   */
float CDecay::getEk3bodyNewton(CFrame* frag1, CFrame* frag2, CFrame* frag3)
{
  float totMass = frag1->A + frag2->A + frag3->A;
  
  CFrame frag1CM(frag1->A);
  CFrame frag2CM(frag2->A);
  CFrame frag3CM(frag3->A);

  float vcm[3] = {0.};
  float ek = 0;
  for (int i=0;i<3;i++)
    {
      vcm[i] += frag1->v[i]*frag1->A;
      vcm[i] += frag2->v[i]*frag2->A;
      vcm[i] += frag3->v[i]*frag3->A;
      vcm[i] /= totMass;

      frag1CM.v[i] = frag1->v[i] - vcm[i];
      frag2CM.v[i] = frag2->v[i] - vcm[i];
      frag3CM.v[i] = frag3->v[i] - vcm[i];
      
    }
      ek += frag1CM.getEnergy();
      ek += frag2CM.getEnergy();
      ek += frag3CM.getEnergy();

  return ek;
}
//*********************************************************
  /**
   * returns the relative kinetic energy of three fragments in their
   * center-of-mass frame. Relativistic version
   \param frag1 is a pointer to fragment 1 
   \param frag2 is a pointer to fragment 2 
   \param frag3 is a pointer to fragment 3 
   */
float CDecay::getEk3bodyRel(CFrame* frag1, CFrame* frag2, CFrame* frag3)
{
    float totMass = frag1->A + frag2->A + frag3->A;
    CFrame frame(totMass);

  for (int i=0;i<3;i++) 
   {
     frame.pc[i] = frag1->pc[i] + frag2->pc[i] + frag3->pc[i];
   }
  frame.totEnergy = frag1->totEnergy + frag2->totEnergy + frag3->totEnergy;
  frame.getVelocityFromMom();

  for (int i=0;i<3;i++)
    {
    frag6Be->pc[i] = frame.pc[i];
    frag6Be->v[i] = frame.pc[i]/frame.pcTot*frame.velocity;
    }
  frag6Be->totEnergy = frame.totEnergy;

  CFrame p1(frag1->A);
  CFrame p2(frag2->A);
  CFrame p3(frag3->A);
  for (int i=0;i<3;i++)
    {
      p1.v[i] = frag1->v[i];
      p2.v[i] = frag2->v[i];
      p3.v[i] = frag3->v[i];
   
    }

  p1.transformVelocity(frame.v);
  p2.transformVelocity(frame.v);
  p3.transformVelocity(frame.v);
  

  ErelRecon = p1.getEnergy() + p2.getEnergy() + p3.getEnergy();

  return ErelRecon;

}





//*******************
float CDecay::getEk3()
{
  return getEk3bodyRel(recon[0],recon[1],recon[2]);
}
//************************************************************


void CDecay::getJacobi(CFrame**part,bool com)
{
  if(Nfrags!=3)
    {
      cout << "what are you doing trying to calculate a jacobi coordinate without 3-particle decay?" << endl;
      return;
    }

  float totMass = part[0]->A + part[1]->A + part[2]->A;

  
  CFrame p1(part[0]->A);
  CFrame p2(part[1]->A);
  CFrame p3(part[2]->A);
  float Etot;
  for (int i=0;i<Nfrags;i++)
    {
      p1.v[i] = part[0]->v[i];
      p2.v[i] = part[1]->v[i];
      p3.v[i] = part[2]->v[i];
   
    }

    if (com)
      {
       CFrame frame(totMass);


      for (int i=0;i<Nfrags;i++) 
       {
         frame.pc[i] = part[0]->pc[i] + part[1]->pc[i] + part[2]->pc[i];
       }
      frame.totEnergy = part[0]->totEnergy + part[1]->totEnergy + 
      part[2]->totEnergy;
      frame.getVelocityFromMom();



      p1.transformVelocity(frame.v);
      p2.transformVelocity(frame.v);
      p3.transformVelocity(frame.v);

      Etot = p1.getEnergy() + p2.getEnergy() + p3.getEnergy();
      }
    else
      {
	Etot = 0.;
	for (int i=0;i<3;i++) Etot += 0.5*p1.A*pow(p1.v[i]/.9784,2)
	  + 0.5*p2.A*pow(p2.v[i]/.9784,2) + 0.5*p3.A*pow(p3.v[i]/.9784,2);
      }
  
  //p-p relative velocity can use newton
  float vrel[3];
  float vrelative = 0.;
  float dot = 0.;
  float v3 = 0.;
  for (int i=0;i<3;i++)
    {
      vrel[i] = p1.v[i] - p2.v[i];
      vrelative += pow(vrel[i],2);
      dot += vrel[i]*p3.v[i];
      v3 += pow(p3.v[i],2);
    }
  vrelative = sqrt(vrelative);
  v3 = sqrt(v3);

  float Epp = 0.5*0.5*pow(vrelative/.9784,2);
  x_T = Epp/Etot;

  CosTheta_T = dot/v3/vrelative;


  //6Be-p1
  vrelative = 0.;
  dot = 0.;
  float v2 = 0.;
  for (int i=0;i<3;i++)
    {
      vrel[i] = p1.v[i] - p3.v[i];
      vrelative += pow(vrel[i],2);
      dot += vrel[i]*p2.v[i];
      v2 += pow(p2.v[i],2);
    }
  vrelative = sqrt(vrelative);
  v2 = sqrt(v2);

  float Epa = 0.5*28./29.*pow(vrelative/.9784,2);
  x_Y[0] = Epa/Etot;

  CosTheta_Y[0] = dot/v2/vrelative;


  //6Be-p2
  vrelative = 0.;
  dot = 0.;
  float v1 = 0.;
  for (int i=0;i<3;i++)
    {
      vrel[i] = p2.v[i] - p3.v[i];
      vrelative += pow(vrel[i],2);
      dot += vrel[i]*p1.v[i];
      v1 += pow(p1.v[i],2);
    }
  vrelative = sqrt(vrelative);
  v1 = sqrt(v1);

  Epa = 0.5*28./29.*pow(vrelative/.9784,2);
  x_Y[1] = Epa/Etot;

  CosTheta_Y[1] = dot/v1/vrelative;
}


void CDecay::getJacobiSecondary()
{
  getJacobi(recon,(bool)1);
}

void CDecay::getJacobiPrimary()
{
  getJacobi(real,(bool)0);
} 



