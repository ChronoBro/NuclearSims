#include "plf.h"
#include <iostream>
#include "TH1S.h"
#include "TH1I.h"
#include "TH2S.h"
#include "TFile.h"
#include "frag.h"
#include "decay.h"
#include "TRandom.h"
#include "MathGeoLib.h"
#include "detectorV2.h"
#include <vector>

using namespace std;


int main()
{
  int correct_ID = 0;
  int total_ID = 0;

  //ofstream ofile("lifetimes.txt");

  bool einstein = 1; //switch for newtonian(0) or relativistic(1) 

  string Loss1("../loss/Hydrogen.loss");
  string Loss2("../loss/Hydrogen.loss");
  string Loss3("../loss/Be/Argon.loss");
  string Loss4("../loss/Be/Calcium.loss");

  
  float EPA = 100.;
  int plfMass = 26;
  float density = 1848.; //mg/cm^3 of berrylium
  
  int const Nevents = 3000000;
  //int const Nevents =  1000;

  //float thickness = 188; //mg/cm2
  //float thickness = 94;
  float thickness = 396;
  //float thickness = 0.0001;  
  float const targetSize = 10.;
  float Csi_Res = .1275;
  //float straggle = 1.01;
  float straggle = 0.;
  //float CsiRes = .113;

  CPlf plf(thickness);
  CFrag frag1(1.,1.,Loss1,Csi_Res,thickness);
  CFrag frag2(1.,1.,Loss2,Csi_Res,thickness);
  CFrag frag3(8.,24.,Loss3,Csi_Res,thickness);

  float zdist=100.;


  float offset = 2*2.5*2.54;

  int nDet = 64;
  vector<detectorV2> neutronArray;

  //figure our scheme for tagging neutrons and correlation between separation of detectors
  for(int i=0;i<4;i++)
    {
      for(int j=0;j<4;j++)
	{
	  neutronArray.push_back(detectorV2(i*2.5*2.54 - offset,j*2.5*2.54 - offset,zdist,12*2.54));	  
	}
    }

   for(int i=0;i<4;i++)
    {
      for(int j=0;j<4;j++)
	{
	  neutronArray.push_back(detectorV2( i*2.5*2.54+1.25*2.54 - offset, j*2.5*2.54 + 1.25*2.54 - offset,zdist,12*2.54));	  
	}
     }

   for(int i=0;i<4;i++)
    {
      for(int j=0;j<4;j++)
	{
	  neutronArray.push_back(detectorV2( i*2.5*2.54+1.25*2.54 - offset, j*2.5*2.54 - offset,zdist + 12*2.54,12*2.54));	  
	}
    }

   for(int i=0;i<4;i++)
     {
       for(int j=0;j<4;j++)
	 {
	   neutronArray.push_back(detectorV2( i*2.5*2.54 - offset, j*2.5*2.54 + 1.25*2.54 - offset,zdist + 12*2.54,12*2.54));	  
	 }
     }

  
  for(int i=0;i<neutronArray.size();i++)
    {
      neutronArray[i].setGeometry(1,1,1,0,1,0,1);
    }


  
  
  // neutronArray.push_back(detectorV2(0.5*2.54,0.5*2.54,zdist,12*2.54));
  // neutronArray.push_back(detectorV2(-0.5*2.54,0.5*2.54,zdist,12*2.54));
  // neutronArray.push_back(detectorV2(0.5*2.54,-0.5*2.54,zdist,12*2.54));
  // neutronArray.push_back(detectorV2(-0.5*2.54,-0.5*2.54,zdist,12*2.54));

  // neutronArray.push_back(detectorV2(1.5*2.54,1.5*2.54,zdist,12*2.54));
  // neutronArray.push_back(detectorV2(-1.5*2.54,1.5*2.54,zdist,12*2.54));
  // neutronArray.push_back(detectorV2(1.5*2.54,-1.5*2.54,zdist,12*2.54));
  // neutronArray.push_back(detectorV2(-1.5*2.54,-1.5*2.54,zdist,12*2.54));

  // neutronArray.push_back(detectorV2(1.5*2.54,0.5*2.54,zdist,12*2.54));
  // neutronArray.push_back(detectorV2(-1.5*2.54,0.5*2.54,zdist,12*2.54));
  // neutronArray.push_back(detectorV2(1.5*2.54,-0.5*2.54,zdist,12*2.54));
  // neutronArray.push_back(detectorV2(-1.5*2.54,-0.5*2.54,zdist,12*2.54));
  
  // neutronArray.push_back(detectorV2(0.5*2.54,1.5*2.54,zdist,12*2.54));
  // neutronArray.push_back(detectorV2(-0.5*2.54,1.5*2.54,zdist,12*2.54));
  // neutronArray.push_back(detectorV2(0.5*2.54,-1.5*2.54,zdist,12*2.54));
  // neutronArray.push_back(detectorV2(-0.5*2.54,-1.5*2.54,zdist,12*2.54));

  
  // neutronArray[0].setGeometry(1,1,1,0,1,0,1);
  // neutronArray[1].setGeometry(1,1,1,0,1,0,1);
  // neutronArray[2].setGeometry(1,1,1,0,1,0,1);
  // neutronArray[3].setGeometry(1,1,1,0,1,0,1);
  // neutronArray[4].setGeometry(1,1,1,0,1,0,1);
  // neutronArray[5].setGeometry(1,1,1,0,1,0,1);
  // neutronArray[6].setGeometry(1,1,1,0,1,0,1);
  // neutronArray[7].setGeometry(1,1,1,0,1,0,1);

  // neutronArray[8].setGeometry(1,1,1,0,1,0,1);
  // neutronArray[9].setGeometry(1,1,1,0,1,0,1);
  // neutronArray[10].setGeometry(1,1,1,0,1,0,1);
  // neutronArray[11].setGeometry(1,1,1,0,1,0,1);
  // neutronArray[12].setGeometry(1,1,1,0,1,0,1);
  // neutronArray[13].setGeometry(1,1,1,0,1,0,1);
  // neutronArray[14].setGeometry(1,1,1,0,1,0,1);
  // neutronArray[15].setGeometry(1,1,1,0,1,0,1);

  
  
  
  detector S4(zdist,0,0); //distance from target along z-axis and thickness (if thickness=0 it speeds things up), and coordinate(not functioning now
  S4.setGeometry(128,128,7.5,62.5,0,0,0); //zeros at end signify it is a ring

  detector CsI1(zdist+5.,39.,0);
  CsI1.setGeometry(1,16,30.5,75,0,0,0);
  //CsI1.setGeometryBack(1,8,42,62.5);

  detector CsI2(zdist+5.,39.,0.);
  CsI2.setGeometry(1,8,7.5,30,0,0,0);
  //CsI2.setGeometryBack(1,4,7.5,42);
  
  CDecay decay(&frag1,&frag2,&frag3,einstein);

  float Qbreakup = decay.EkTot8B;

  TH1I JacobiT_x_p("JacobiT_x_p","",10,0,1);
  TH1I JacobiY_x_p("JacobiY_x_p","",10,0,1);
  TH1I JacobiT_y_p("JacobiT_y_p","",10,-1,1);
  TH1I JacobiY_y_p("JacobiY_y_p","",10,-1,1);
  
  TH1I JacobiT_x_s("JacobiT_x_s","",10,0,1);
  TH1I JacobiY_x_s("JacobiY_x_s","",10,0,1);
  TH1I JacobiT_y_s("JacobiT_y_s","",20,-1,1);
  TH1I JacobiY_y_s("JacobiY_y_s","",20,-1,1);
  TH2S JacobiT_xy_s("JacobiT_xy_s","",10,0,1,20,-1,1);
  TH2S JacobiY_xy_s("JacobiY_xy_s","",10,0,1,20,-1,1);
  TH1I JacobiT_y1_s("JacobiT_y1_s","",20,-1,1);
  TH1I JacobiY_y1_s("JacobiY_y1_s","",20,-1,1);
  TH1I JacobiT_y2_s("JacobiT_y2_s","",20,-1,1);
  TH1I JacobiY_y2_s("JacobiY_y2_s","",20,-1,1);
  TH1I JacobiT_y3_s("JacobiT_y3_s","",20,-1,1);
  TH1I JacobiY_y3_s("JacobiY_y3_s","",20,-1,1);
  
  TH1I hist_vel("vel","vel",100,0.,14.);
  TH1I hist_theta("theta","theta",200,0,10);
  TH1I hist_phi("phi","phi",125,0,360);

  TH1I hist_vel_R("velR","vel",100,0.,14.);
  TH1I hist_theta_R("thetaR","theta",200,0,10);
  TH1I hist_phi_R("phiR","phi",125,0,360);
  TH1I hist_theta_emission_R("thetaEmissionR","thetaEmission",100,0,180);

  TH1I alpha_mom("alpha_mom","",200,0,200);

  TH1I keAlpha("keAlpha","",1000,0,1000);
  TH1I keProton("keProton","",1000,0,200);

  TH1I hist_aRatio("aRatio","",200,0,1);

  TH1S hist_Erel("Erel","Erel",100,0.,10.);
  TH1S Erel_noT("Erel_noT","Erel_noT",500,0.,10.);
  TH1S hist_ppRel("ppRel","",500,-3,10);
  TH2S map_Erel_dthick("Erel_dthick","Erel_dthick",200,1.5,4.5,100,0,1000);
  TH2S map_Erel6Be_dthick("Erel6Be_dthick","Erel6Be_dthick",200,-1.,1.,100,0,1000);

  TH2S protonXY_s("protonXY_s","protonxy",250,-20,20,250,-13,13);
  TH2S protonXY_s_1("protonXY_s_1","protonxy",250,-20,20,250,-13,13);
  TH2S protonXY_s_2("protonXY_s_2","protonxy",250,-20,20,250,-13,13);
  TH2S protonXY_s_1o("protonXY_s_1o","protonxy",250,-20,20,250,-13,13);
  TH2S protonXY_s_2o("protonXY_s_2o","protonxy",250,-20,20,250,-13,13);
  TH2S protonXY_s_1i("protonXY_s_1i","protonxy",250,-20,20,250,-13,13);
  TH2S protonXY_s_2i("protonXY_s_2i","protonxy",250,-20,20,250,-13,13);
  TH2I protonCsI("protonCsI","protonCsI",16,-0.5,15.5,16,-0.5,15.5);

  TH2S coreXY_s("coreXY_s","alphaxy",250,-20,20,250,-13,13);

  TH1I hist_momDisttr("momDisttr","",100,-400,400);
  TH1I hist_momDistz("momDistz","",100,-400,400);

  TH1F ArEdist("ArEdist","",1000,10000,12000);
  
  TRandom lifetime_distribution(0);


  float mean =0.;
  float sig = 0.;

  double decayconstant =15.0; //in ns
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

  CLoss CaInBe(Loss4,35);
  CLoss ArInBe(Loss3,32);
  float ECa = EPA*35;

  float3 origin = origin.zero;

  TRandom tStart;
  
  for (int i=0;i<Nevents;i++)
    {

      float EPA0=0;
      EPA0=EPA;

      
      if (i%1000 == 0) cout << '\xd' <<  i << flush;
      dthick = thickness*plf.ran.Rndm();
      plf.GetPlf(EPA0,26,1);
      // double EoutCa = CaInBe.getEout(ECa,dthick);

      // EPA0 = EoutCa/(float)35;
      // plf.MultiScat(dthick);
      //cout << "EPA of 35Ca = " << EPA0 << endl;

      xTarget = (1.-2.*plf.ran.Rndm())/2.*targetSize;
      yTarget = (1.-2.*plf.ran.Rndm())/2.*targetSize;
      plf.setBeamSpot(targetSize);
      
      plf.GetPlf(EPA0,24,0);

      //lifetime = (float)lifetime_distribution.Exp(plf.frame->gamma*decayconstant);
      // xTarget = xTarget + plf.frame->v[0]*lifetime*10.; //multiplied by 10 to put in mm
      // yTarget = yTarget + plf.frame->v[1]*lifetime*10.;
      // zBreakup = lifetime*plf.frame->v[2]*10. + dthick/density*10.;

      zBreakup = dthick/density*10.;

      //plf.GetPlf();



 
      double EoutCa2;
     

 

 
      decay.ModeMicroCanonical();
      
      frag1.AddVelocity(plf.frame->v); //need to transform these right after decay so that the energies in the LAB frame are correct
      frag2.AddVelocity(plf.frame->v);
      frag3.AddVelocity(plf.frame->v);

 

      
      hist_vel.Fill(plf.frame->velocity);
      hist_theta.Fill(plf.frame->theta*180./plf.pi);
      hist_phi.Fill(plf.frame->phi*180./plf.pi);
      frag3.real->getVelocity();
      //cout << frag3.real->pcTot << endl;
      ArEdist.Fill(frag3.real->pcTot);      
      
      //cout << "thickness in mg/cm2 = " << thickness << endl;
      //cout << "Zbreakup in mg/cm2 = " << zBreakup/10.*density << endl;
      //cout << "frag3 energy = " << frag3.real->energy << endl;
      //EPA0 = frag3.real->energy/32.;
      
     
      
      //multiple scattering of plf
      //if (thickness > 0.)  plf.MultiScat(thickness-dthick);

      

      JacobiT_x_p.Fill(decay.x_T);
      JacobiT_y_p.Fill(decay.CosTheta_T);
      JacobiT_y_p.Fill(-decay.CosTheta_T);
 
      JacobiY_x_p.Fill(decay.x_Y[0]);
      JacobiY_x_p.Fill(decay.x_Y[1]);
      JacobiY_y_p.Fill(-decay.CosTheta_Y[0]);
      JacobiY_y_p.Fill(-decay.CosTheta_Y[1]);


 
      frag1.setTrajectory(xTarget/10., yTarget/10., zBreakup/10.);
      frag2.setTrajectory(xTarget/10., yTarget/10., zBreakup/10.);
      int l=0;
      int m=0;
      int nhit=frag3.hit4(xTarget,yTarget,zBreakup);
      
      
      float start_res = 1.; //ns
      float startTime = tStart.Gaus(0,start_res/2.355);

      
      //will ultimately need some logic to decide on interaction if neutron trajectory is through 2 (or more) detectors
      for( l=0;l<neutronArray.size();l++)
	{
	  if( neutronArray[l].event(frag1.trajectory) )
	    {
	      nhit++;
	      for( m=0;m<neutronArray.size();m++)
		{
		  if( neutronArray[m].event(frag2.trajectory) )
		    {
		      frag1.setRecon(neutronArray[l].distInteraction, neutronArray[l].measuredDist, startTime);
		      frag2.setRecon(neutronArray[m].distInteraction, neutronArray[m].measuredDist, startTime);
		      nhit++;
		      break;
		    }
		}
	    }
	}

 
      if (nhit != 3 || (l==m))
	{
	  continue;
	}
      
      Ndet++;

      keProton.Fill(frag1.real->energy);
      keProton.Fill(frag2.real->energy);
      keProton.Fill(frag3.real->energy);

 
      NdetClean++;

      //correct for energy loss in half of target and get velocity

      float ErelR = decay.getErelRecon();
      Erel_noT.Fill(ErelR);
      
      // frag1.Egain(thickness/2.);
      // frag2.Egain(thickness/2.);
      frag3.Egain(thickness/2.);

      ErelR = decay.getErelRecon();
      hist_Erel.Fill(ErelR);     


      //float ErelR = decay.getEk5();

      hist_theta_R.Fill(decay.plfRecon->theta*180./3.1459);
      hist_vel_R.Fill(decay.plfRecon->velocity);




      for (int ii=0;ii<6;ii++)hist_ppRel.Fill(decay.ppRel[ii]);
      hist_aRatio.Fill(decay.aRatio);

      map_Erel_dthick.Fill(ErelR,dthick);
      hist_theta_R.Fill(180. - decay.plfRecon->theta*180./plf.pi);

      decay.getJacobiSecondary();

      JacobiT_x_s.Fill(decay.x_T);
      JacobiT_y_s.Fill(decay.CosTheta_T);
      JacobiT_y_s.Fill(-decay.CosTheta_T);
      JacobiT_xy_s.Fill(decay.x_T,decay.CosTheta_T);
      JacobiT_xy_s.Fill(decay.x_T,-decay.CosTheta_T);
      if (decay.x_T < 0.33333)
	{
	  JacobiT_y1_s.Fill(decay.CosTheta_T);
	  JacobiT_y1_s.Fill(-decay.CosTheta_T);
	}
      else if (decay.x_T < 0.666666)
	{
	  JacobiT_y2_s.Fill(decay.CosTheta_T);
	  JacobiT_y2_s.Fill(-decay.CosTheta_T);
	}
      else 
	{
	  JacobiT_y3_s.Fill(decay.CosTheta_T);
	  JacobiT_y3_s.Fill(-decay.CosTheta_T);
	}

      //cout << "decay.Etot: " << decay.Etot << endl;
      //cout << "checking frag1->energy: " << decay.Ep << endl;
      //cout << "checking Jacobi Y: " << decay.frag[0]->energy/decay.Etot << endl;
      JacobiY_x_s.Fill(decay.Ep/decay.Etot);
      JacobiY_x_s.Fill(decay.Ep2/decay.Etot);
      JacobiY_y_s.Fill(-decay.CosTheta_Y[0]);
      JacobiY_y_s.Fill(-decay.CosTheta_Y[1]);
      JacobiY_xy_s.Fill(decay.x_Y[0],-decay.CosTheta_Y[0]);
      JacobiY_xy_s.Fill(decay.x_Y[1],-decay.CosTheta_Y[1]);

      if (decay.x_Y[0] < 0.33333) JacobiY_y1_s.Fill(-decay.CosTheta_Y[0]);
      else if (decay.x_Y[0] < 0.6666) JacobiY_y2_s.Fill(-decay.CosTheta_Y[0]);
      else JacobiY_y3_s.Fill(-decay.CosTheta_Y[0]);


      if (decay.x_Y[1] < 0.33333) JacobiY_y1_s.Fill(-decay.CosTheta_Y[1]);
      else if (decay.x_Y[1] < 0.6666) JacobiY_y2_s.Fill(-decay.CosTheta_Y[1]);
      else JacobiY_y3_s.Fill(-decay.CosTheta_Y[1]);

      float x = frag1.recon->theta*180./plf.pi*cos(frag1.recon->phi);
      float y = frag1.recon->theta*180./plf.pi*sin(frag1.recon->phi);
      protonXY_s.Fill(x,y);
      protonXY_s_1.Fill(x,y);
      protonCsI.Fill(frag1.Ring->CsIstrip,frag2.Ring->CsIstrip);
      if (frag1.Ring->CsIring == 0)
	{
	  protonXY_s_1i.Fill(x,y);
	}
      else
	{
	  protonXY_s_1o.Fill(x,y);
	}

      x = frag2.recon->theta*180./plf.pi*cos(frag2.recon->phi);
      y = frag2.recon->theta*180./plf.pi*sin(frag2.recon->phi);
      protonXY_s.Fill(x,y);
      protonXY_s_2.Fill(x,y);
      if (frag2.Ring->CsIring == 0)
	{
	  protonXY_s_2i.Fill(x,y);
	}
      else
	{
	  protonXY_s_2o.Fill(x,y);
	}


      x = frag3.recon->theta*180./plf.pi*cos(frag3.recon->phi);
      y = frag3.recon->theta*180./plf.pi*sin(frag3.recon->phi);
      coreXY_s.Fill(x,y);

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

  TFile f ("../root/simO26.root","RECREATE");


  hist_vel.Write();
  hist_theta.Write();
  hist_phi.Write();
  hist_Erel.Write();
  Erel_noT.Write();
  hist_theta_R.Write();
  hist_vel_R.Write();


  keAlpha.Write();
  keProton.Write();
  ArEdist.Write();
  alpha_mom.Write();
  protonXY_s.Write();
  protonXY_s_1.Write();
  protonXY_s_2.Write();
  protonXY_s_1i.Write();
  protonXY_s_1o.Write();
  protonXY_s_2i.Write();
  protonXY_s_2o.Write();
  protonCsI.Write();
  coreXY_s.Write();

  map_Erel_dthick.Write();
  map_Erel6Be_dthick.Write();
  hist_ppRel.Write();
  hist_aRatio.Write();


  JacobiT_x_s.Write();
  JacobiT_y_s.Write();
  JacobiT_y1_s.Write();
  JacobiT_y2_s.Write();
  JacobiT_y3_s.Write();
  JacobiY_x_s.Write();
  JacobiY_y_s.Write();
  JacobiY_y1_s.Write();
  JacobiY_y2_s.Write();
  JacobiY_y3_s.Write();
  JacobiT_xy_s.Write();
  JacobiY_xy_s.Write();

  JacobiT_x_p.Write();
  JacobiT_y_p.Write();
  JacobiY_x_p.Write();
  JacobiY_y_p.Write();



  //ofile.close();

  f.Write();

}
