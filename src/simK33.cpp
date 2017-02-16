#include "plf.h"
#include <iostream>
#include "TH1S.h"
#include "TH1I.h"
#include "TH2S.h"
#include "TFile.h"
#include "frag.h"
#include "decay.h"
#include "TRandom.h"
//#include "CMomDist.h"


using namespace std;


int main()
{
  int correct_ID = 0;
  int total_ID = 0;

  //ofstream ofile("lifetimes.txt");

  bool einstein = 1; //switch for newtonian(0) or relativistic(1) 

  string Loss1("../loss/Hydrogen.loss");
  string Loss2("../loss/Hydrogen.loss");
  string Loss3("../loss/Argon.loss");
  

  float plfMass = 33.;
  float targetMass = 9.;
  
  int const Nevents = 50000;
  //int const Nevents =  1000;

  //float thickness = 188; //mg/cm2
  //float thickness = 94;
  float thickness = 396;
  float const targetSize = 10.;
  float Csi_Res = .1275;
  //float straggle = 1.01;
  float straggle = 0.;
  //float CsiRes = .113;

  float zdist=300.;
  
  float EPA0 = 70;  
  //energy loss in 1/2 of target
  //float epa2_mean = Ploss.getEout(EPA0*plfMass,thickness/2.)/plfMass;
  float Pbeam2_mean = sqrt(pow((EPA0+931.478)*plfMass,2)-pow(plfMass*931.478,2));


  //reaction center of mass velocity
  float Ecm = (EPA0+931.478)*plfMass + targetMass*931.478;
  float pcm = Pbeam2_mean;
  float vCM = pcm/Ecm*30.;

  


  CPlf plf(thickness);
  //CFrag frag1(1.,1.,Loss1,Csi_Res,thickness);
  CFrag frag2(1.,1.,Loss2,Csi_Res,thickness);
  CFrag frag3(18.,32.,Loss3,Csi_Res,thickness);

 

  detector S4(zdist,0,0); //distance from target along z-axis and thickness (if thickness=0 it speeds things up), and coordinate(not functioning now
  S4.setGeometry(128,128,7.5,62.5,0,0,0); //zeros at end signify it is a ring
  detector CsI1(zdist+5.,39.,0);
  CsI1.setGeometry(1,8,35,62.5,0,0,0);
  //CsI1.setGeometryBack(1,8,42,62.5);
  detector CsI2(zdist+5.,39.,0.);
  CsI2.setGeometry(1,4,7.5,35,0,0,0);
  //CsI2.setGeometryBack(1,4,7.5,42);
  
  CDecay decay(&frag2,&frag3,einstein);

  float Qbreakup = decay.EkTot8B;
  
  TH1I hist_vel("vel","vel",100,0.,14.);
  TH1I hist_theta("theta","theta",200,0,10);
  TH1I hist_phi("phi","phi",125,0,360);

  TH1I hist_vel_R("velR","vel",100,0.,6.);
  TH1I hist_theta_R("thetaR","theta",200,0,10);
  TH1I hist_phi_R("phiR","phi",125,0,360);
  TH1I hist_theta_emission_R("thetaEmissionR","thetaEmission",100,0,180);

  TH1I alpha_mom("alpha_mom","",200,0,200);

  TH1I keAlpha("keAlpha","",1000,0,1000);
  TH1I keProton("keProton","",1000,0,200);

  TH1I hist_aRatio("aRatio","",200,0,1);

  TH1S hist_Erel("Erel","Erel",600,0.,10.);
  TH1S hist_ppRel("ppRel","",500,-3,10);
  TH2S map_Erel_dthick("Erel_dthick","Erel_dthick",200,1.5,4.5,100,0,1000);
  TH2S map_Erel6Be_dthick("Erel6Be_dthick","Erel6Be_dthick",200,-1.,1.,100,0,1000);

  TH2F hist_Ep_cosTheta("Ep_cosTheta","Ep_cosTheta",100,-1,1,50,0,5);

  TH2S coreXY_s("coreXY_s","alphaxy",250,-20,20,250,-13,13);

  TH1I hist_momDisttr("momDisttr","",100,-400,400);
  TH1I hist_momDistz("momDistz","",100,-400,400);

  TRandom lifetime_distribution(0);


  float mean =0.;
  float sig = 0.;

  double decayconstant = 0.1; //in ns
  float lifetime = 0;

  int Ndet = 0;
  int NdetClean = 0;
  int NdetCleaner = 0;
  int NdetCleaner6 = 0;
  int NdetCleanerCen = 0;
  int Nmissingp = 0;
  int nguardring = 0;
  int highenergyprotons = 0;

  float dthick = 0.;
  float xTarget = 0.;
  float yTarget = 0.;
  float zBreakup = 0.;

  cout << endl << "Pbeam set to: " << Pbeam2_mean << endl << endl;

  for (int i=0;i<Nevents;i++)
    {

      if (i%1000 == 0) cout << '\xd' <<  i << flush;
      dthick = thickness*plf.ran.Rndm();
      xTarget = (1.-2.*plf.ran.Rndm())/2.*targetSize;
      yTarget = (1.-2.*plf.ran.Rndm())/2.*targetSize;

   
      plf.GetPlf(Pbeam2_mean,plfMass);
      //multiple scattering of plf
      if (thickness > 0.)  plf.MultiScat(1.-dthick/thickness);

      hist_vel.Fill(plf.frame->velocity);
      hist_theta.Fill(plf.frame->theta*180./plf.pi);
      hist_phi.Fill(plf.frame->phi*180./plf.pi);

      plf.frame->getVelocity();
      
      //lifetime = (float)lifetime_distribution.Exp(plf.frame->gamma*decayconstant);
      //lifetime = 0.;

      
      
      xTarget = xTarget + plf.frame->v[0]*lifetime*10.; //divided by 10 to put in mm
      yTarget = yTarget + plf.frame->v[1]*lifetime*10.;
      zBreakup = lifetime*plf.frame->v[2]*10.;

      //ofile << lifetime << " " << zBreakup << endl;

      //now break up plf
      plf.isotropic();

      //cout << endl << "before decay mode frag2.real->theta = " << frag2.real->theta << endl;
      
      //decay.ModeMoscow(); //sequential micro+P26Be
      //decay.ModeMicroCanonical();
      decay.Mode(1.96,0.,1.,32.);

       //frag1.AddVelocity(plf.frame->v);
      frag2.AddVelocity(plf.frame->v);
      frag3.AddVelocity(plf.frame->v);

      //cout << endl << "after decay mode frag2.real->theta = " << frag2.real->theta << endl;      


      //interaction in target, continue if stopped
      //if (frag1.targetInteraction(dthick,thickness)) continue;
      if (frag2.targetInteraction(dthick,thickness)) continue;
      if (frag3.targetInteraction(dthick,thickness)) continue;

      // cout << endl<< "original frag1 theta = " << frag1.real->theta << endl;
      // cout << "original frag2 theta = " << frag2.real->theta << endl;
      // cout << "original frag1 energy = " << frag1.real->energy << endl;
      // cout << "original frag2 energy = " << frag2.real->energy << endl;
      
      //int nhit = frag1.hit3(xTarget, yTarget, zBreakup, straggle) + frag2.hit3(xTarget, yTarget, zBreakup,straggle);
      //nhit += frag3.hit4(xTarget, yTarget, zBreakup);

      //if(nhit!=3) continue;

      int nhit;
      int nhit1 = S4.event(frag2.real->theta, frag2.real->phi,0);//(frag1.real->theta, frag1.real->phi,0);
           
      // int nhit2 = S4.event(frag2.theta_prime, frag3.phi_prime,0);
      // frag2.recon->theta = S4.thetaHit;
      // frag2.recon->phi = S4.phiHit;

      if(nhit1)
      	{
	  frag2.recon->theta = S4.thetaHit;
	  frag2.recon->phi = S4.phiHit;
      	  frag2.recon->energy = frag2.real->energy  + sqrt(frag2.real->energy)*frag2.CsI_res*
      	    frag2.ran.Gaus(0.,1.);
      	  frag2.recon->getVelocity();
	  //cout << "should be similar frag1 energy = " << frag1.recon->energy << endl;
	  //cout << endl << "frag2.real->theta = " << frag2.real->theta << endl;
	  //cout << "frag2.recon->theta = " << frag2.recon->theta << endl;

	}

      
      
      
      // if(nhit2)
      // 	{
      // 	  frag2.recon->energy = frag2.real->energy  + sqrt(frag2.real->energy)*frag2.CsI_res*
      // 	    frag2.ran.Gaus(0.,1.);
      // 	  frag2.recon->getVelocity();
      // 	  //cout << "should be similar frag2 energy = " << frag2.recon->energy << endl;
      // 	}
      
      
      nhit = nhit1 + frag3.hit4(xTarget, yTarget, zBreakup);
      
      
      //if(frag3.real->theta > frag1.Ring->theta_min)
      //if (nhit != 3 || zdist < zBreakup) continue;
      if(nhit!=2) continue;
      
      Ndet++;

      //keProton.Fill(frag1.real->energy);
      keProton.Fill(frag2.real->energy);
      keProton.Fill(frag3.real->energy);

      int checkSeg1;
      int checkSeg2;
      
      if(CsI1.event(frag2.real->theta,frag2.real->phi,0))
      	{
      	  checkSeg1 = CsI1.segmentYhit;
      	  if(CsI1.event(frag3.real->theta, frag3.real->phi,0))
      	    {
      	      checkSeg2 = CsI1.segmentYhit;
      	      if(checkSeg1==checkSeg2)
      		{
      		  continue;
      		}
      	    }
      	  // else if(!CsI2.event(frag3.real->theta, frag3.real->phi,0))
      	  //   {
      	  //     continue;
      	  //   }
      	}
      else if(CsI2.event(frag2.real->theta,frag2.real->phi,0))
      	{
      	  checkSeg1 = CsI2.segmentYhit;
      	  if(CsI2.event(frag3.real->theta,frag3.real->phi,0))
      	    {
      	      checkSeg2 = CsI2.segmentYhit;
      	      if(checkSeg1==checkSeg2)
      		{
      		  continue;
      		}
      	    }
      	  // else if(!CsI1.event(frag3.real->theta, frag3.real->phi,0))
      	  //   {
      	  //     continue;
      	  //   }
      	}
      else
      	{
      	  continue;
      	}

      

	  
	  //if (decay.OnTopOf3()) continue;

      NdetClean++;

      //correct for energy loss in half of target and get velocity
      
      //frag2.Egain(thickness/2.);
      frag2.Egain(thickness/2.);
      frag3.Egain(thickness/2.);


      float ErelR = decay.getErelRecon();
      //cout << "WTF is with ErelR? " << ErelR << endl;

      //float ErelR = decay.getEk5();
      hist_Erel.Fill(ErelR);
      
      hist_theta_R.Fill(decay.plfRecon->theta*180./3.1459);
      hist_vel_R.Fill(decay.plfRecon->velocity);
      hist_Ep_cosTheta.Fill(decay.partCM[0]->v[2]/decay.partCM[0]->velocity,ErelR);
      //cout << "decay.plfRecon->theta*180./3.14159 = " << decay.plfRecon->theta*180./3.1415927 << endl;
      
    
    }

  cout << endl;
  cout << "number with 1 lost proton " << (float)Nmissingp/(float)Nevents << endl;

  cout << "detection efficiency = " << (float)Ndet/(float)Nevents << endl;
  cout << "clean detection efficiency (separate stips and Csi) = " 
       << (float)NdetClean/(float)Nevents << endl;
//   cout << "final detection efficiency (protons stay in Csi) = " 
//        << (float)NdetCleaner/(float)Nevents << endl;
//   cout << "final detection efficiency (protons stat in Csi)6Li in det6 = " 
//        << (float)NdetCleaner6/(float)Nevents << endl;

//   cout << "final detection efficiency (protons stat in Csi)(3towers) = " 
//        << (float)NdetCleanerCen/(float)Nevents << endl;

  cout << "fraction of cores hitting guard ring = "  
       << (float)nguardring/(float)Nevents << endl;
  cout << "fraction of high energy protons being thrown out = "
       << (float)highenergyprotons/(float)Nevents << endl;

  mean = mean/(float)Ndet;
  sig = sig/(float)(Ndet-1) - (float)(Ndet)/(float)(Ndet-1)*pow(mean,2);
  sig = sqrt(sig);
  cout << "mean Erel= " << mean << " mean Ex " << mean+Qbreakup 
       << " sig= " << sig << " FWHM= " 
       << sig*2.35 << endl;
  

  cout << "correct identification of 6B3 p-p is " << 
    (float)correct_ID/(float)total_ID << " events= " <<  
    (float)total_ID/(float)NdetCleaner << endl;

  TFile f ("../root/simK33.root","RECREATE");


  hist_vel.Write();
  hist_theta.Write();
  hist_phi.Write();
  hist_Erel.Write();
  hist_theta_R.Write();
  hist_vel_R.Write();
  hist_Ep_cosTheta.Write();

  keAlpha.Write();
  keProton.Write();
  alpha_mom.Write();

  map_Erel_dthick.Write();
  map_Erel6Be_dthick.Write();
  hist_ppRel.Write();
  hist_aRatio.Write();

  //ofile.close();

  f.Write();

}
