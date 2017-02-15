#include "array.h"

using namespace std;


/**
 * Constructor
 \param dist0 is radius to si detector from  nominal target center. (cm)
\param fromCenter0 offset of real target from nominal target center (cm)
 \param Active is the active dimension of Si detector (cm)
 \param Physical is the physical dimension of Si detector (cm)
 */


CArray::CArray(int Ndetectors0,char * nameMapFile)
{
  Ndetectors = Ndetectors0;
  cout << endl << "Initializing detector array..." << endl;
  mapFile.open(nameMapFile);
  map = new int*[2*Ndetectors];
  
}

int CArray::setMap()
{
  int imax=0;
  int jmax=0;
  for(int i; i < Ndetectors; i++)
    {
      int check = Detectors[i].didSetGeo;
      if(check)
	{
	  imax = max(Detectors[i].Nx,imax);
	  jmax = max(Detectors[i].Ny, jmax);
	  continue;
	}
      else
	{
	  cout << endl << "Didn't set geometry for detector " << i << endl;
	  return 0;
	}
    }

  total = (imax+1)*(jmax+1);
  
  for(int i=0; i < 2*Ndetectors; i++)
    {
      map[i] = new int[total];
    }
  
  if(mapFile.is_open())
    {
      for(int i; i < total; i++)
	{
	  for(int d; d< 2*Ndetectors; d++)
	    {
	      mapFile >> map[d][i]; //data stored as Nxisegment0, Nyisegment0, Nxisegment1, Nyisegment1 etc...
	    }
	}
    }
  else
    {
      cout << endl << "Couldn't open mapfile" << endl ;
      cout << endl << "Need Map file if you want easy clean detection" << endl;
      return 0;
    }

  
  return 1;

}


int CArray::cleanDetection(int detector1,int segment1, int detector2, int segment2, int isysegment)
{
  if(isysegment)
    {
      for(int i=0;i<total;i++)
	{
	  if(segment1==map[2*detector1][i])
	    {
	      if(segment2==map[2*detector2][i])
		{
		  return 1;
		}
	      return 0;
	      break;
	    }
	}
    }
  else
    {
      for(int i=0;i<total;i++)
	{
	  if(segment1==map[detector1][i])
	    {
	      if(segment2==map[detector2][i])
		{
		  return 1;
		}
	      return 0;
	      break;
	    }
	}

    }
  return 0;
  
}


CArray::CArray(float dist0, float fromCenter0, float Active, float Physical)
{
  pi = acos(-1.);

  dist = dist0;
  fromCenter = fromCenter0;
  xActive = Active; //in-plane active dimension of Si detector (cm)
  yActive = Active; //out-of-plane active dimension of Si detector (cm)

  float const deltaAlpha = 7.5/180.*3.14159; //9.75/180.*3.14159;
  float alpha[5]={-15.25,-7.75,0.,7.75,15.25};
  for (int i=0;i<5;i++) alpha[i] *= 3.14159/180.;
  //for (int i=0;i<5;i++) alpha[i] = (float)(i-2)*deltaAlpha;
  float const deltaBeta = 7.5/180.*3.14159; 
  float const startBeta = 5./180.*3.14159;
  //float const startBeta = deltaBeta;
  nTower = 5;
  float beta[nTower][5];



  //find the gap between the Si
  float spanAngle = atan(Physical/2./dist);  //half the angle spans by the physical Si
  float innerEdge = startBeta - spanAngle; // angle of inner edge of inner Si
  float distEdge = sqrt(pow(dist,2)+pow(Physical/2.,2)); //dist to this edge
  float sepEdge = distEdge*sin(innerEdge)*2.; //sepparation between inner edges

  //cout << "separation between inner edges = " << sepEdge << endl; 
 
  /*
  int nTele[] = {2,4,4,4,2};
  //first tower
  beta[0][0] = -deltaBeta/2.;
  beta[0][1] = deltaBeta/2.;
  //second tower
  beta[1][0] = -1.5*deltaBeta;
  beta[1][1] = -0.5*deltaBeta;
  beta[1][2] = 0.5*deltaBeta;
  beta[1][3] = 1.5*deltaBeta;
  //third tower
  beta[2][0] = -startBeta - deltaBeta;
  beta[2][1] = -startBeta;
  beta[2][2] = startBeta;
  beta[2][3] = startBeta + deltaBeta;
  //forth tower
  beta[3][0] = -1.5*deltaBeta;
  beta[3][1] = -0.5*deltaBeta;
  beta[3][2] = 0.5*deltaBeta;
  beta[3][3] = 1.5*deltaBeta;
  //fifth tower
  beta[4][0] = -deltaBeta/2.;
  beta[4][1] = deltaBeta/2.;
  */
  int nTele[] = {2,3,4,3,2};
  //first tower
  beta[0][0] = -deltaBeta/2.;
  beta[0][1] = deltaBeta/2.;
  //second tower
  beta[1][0] = -deltaBeta;
  beta[1][1] = 0.;
  beta[1][2] = deltaBeta;
  //third tower
  beta[2][0] = -startBeta - deltaBeta;
  beta[2][1] = -startBeta;
  beta[2][2] = startBeta;
  beta[2][3] = startBeta + deltaBeta;
  //forth tower
  beta[3][0] = -deltaBeta;
  beta[3][1] = 0.;
  beta[3][2] = deltaBeta;
  //fifth tower
  beta[4][0] = -deltaBeta/2.;
  beta[4][1] = deltaBeta/2.;


  /*
  int nTele[] = {1,5,4,5,1};
  //first tower
  beta[0][0] = 0.;

  //second tower
  beta[1][0] = -2.*deltaBeta;
  beta[1][1] = -deltaBeta;
  beta[1][2] = 0.;
  beta[1][3] = deltaBeta;
  beta[1][4] = 2.*deltaBeta;
  //third tower
  beta[2][0] = -startBeta - deltaBeta;
  beta[2][1] = -startBeta;
  beta[2][2] = startBeta;
  beta[2][3] = startBeta + deltaBeta;
  //forth tower
  beta[3][0] = -2.*deltaBeta;
  beta[3][1] = -deltaBeta;
  beta[3][2] = 0.;
  beta[3][3] = deltaBeta;
  beta[3][4] = 2.*deltaBeta;
  //fifth tower
  beta[4][0] = 0.;
  */


  float rad_si = 50.;
  float extra_dist = 8.;
  for (int itower=0;itower<nTower;itower++)
    {

     Tower[itower] = new CTower(nTele[itower],dist,fromCenter,alpha[itower],
      &beta[itower][0],xActive,yActive);
    }
}
//**********************************************
  /**
   * destructor
   */
CArray::~CArray()
{
  for (int i=0;i<nTower;i++) delete Tower[i];
  //delete Detectors;
  mapFile.close();
  delete map;
}
//**********************************************
  /**
   * logical indicating if a fragment was detected in the Si array
   \param theta is the zenith angle of fragment in radians
   \param phi is the azimuth angle of particle in radians
   \param xtarget in-plane location of interaction in target from nomimal beam center (cm)
   \param ytarget out-of-plabe location (cm)

   */
int CArray::hit(float theta, float phi, float xtarget, float ytarget)
{
  int ihit =0;
  for (int itower =0;itower<nTower;itower++)
    {
      ihit += Tower[itower]->hit(theta, phi, xtarget, ytarget);
      if (ihit == 1)
	{
	  hitTower = itower;
          thetaRecon = Tower[itower]->thetaRecon;
          phiRecon = Tower[itower]->phiRecon;
          return ihit;
	}
    }
  return ihit;
}
