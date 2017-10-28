#include "detector.h"


/**
 * Constructor
 \param Z0 is the atomic number of fragment
 \param mass0 is mass of fragment in amu
 \param filename specifies the name of the file with energy loss information
 \param CsI_res0 fractional energy resolution of this fragment in CsI
 \param thickness is target thickness in mg/cm2 
*/
using namespace std;

detector::detector(float dist0, float thickness0, int coordinate0)
{
  thickness = thickness0;
  coordinate = coordinate0; //0 to use spherical coordinates, 1 to use cartesian
  didSetGeo=0;
  xD = 0.;
  yD=0.;
  zD=dist0;
  double p0[3] = {0,0,dist0};
  double n0[3] = {0,0,-1};
  detPlane.setPlaneCenter(p0);

  //by default will set plane coordinates to correspind to x-y axes and normal vector along z-axis
  //need to make this more robust
  setDetNormal(n0);
  
  p0[0] = 1.;
  n0[1] = 1.;
  n0[2] = dist0;
  
  detPlane.setPlaneCoordinates(p0,n0);
  
  
}

detector::detector(float x0, float y0, float z0, float thickness0, int coordinate0)
{
  thickness = thickness0;
  coordinate = coordinate0;
  didSetGeo=0;
  xD = x0;
  yD = y0;
  zD = z0;
  double p0[3] = {x0,y0,z0};
  detPlane.setPlaneCenter(p0);
}

int detector::setGeometry(int nRing, int Npie, float rmin, float rmax)
{
  Nx = nRing;
  Ny = Npie;
  
 
  setGeometryFront(nRing,Npie,rmin,rmax);

  if(thickness >0)
    {
      setGeometryBack(nRing,Npie,rmin,rmax);
    }

  
  return 1;
}

int detector::findDetNormal(double* &p1, double* &p2, double* &p3)
{
  detPlane.findNormal(p1,p2,p3);
  return 1;
}


int detector::setDetNormal(double * n)
{
  detPlane.setNormal(n);

  if(thickness > 0)
    {

      double p[3]={xD,yD,zD};
      //adds thickness to plane center
      for(int i=0;i<3;i++)
	{
	  p[i] += n[i]*-1.*thickness;
	}
      
      detBPlane.setNormal(n);
      detBPlane.setPlaneCenter(p);
    }

  return 1;
}

double * detector::labCoordinates()
{
  return detPlane.labCoordinate(xHit,yHit);
}

void detector::moveCenter(double x0, double y0, double z0)
{
  detPlane.moveCenter(x0,y0,z0);
  xD=x0;
  yD=y0;
  zD=z0;
  if(thickness > 0)
    {
      detBPlane.moveCenter(x0,y0,z0);
    }
}

int detector::setGeometryFront(int nRing, int Npie, float rmin, float rmax)
{
  dr = (rmax-rmin)/(float)nRing;
  
  for(int iring=0; iring<nRing; iring++)
    {
      for(int ipie=0;ipie< Npie; ipie++)
	{

	  Location[iring][ipie].fx = xD;
	  Location[iring][ipie].fy =  yD +rmin + dr*(float)(iring+1); //lab coordinates will be set to upper edge of each ring in detector frame
	  Location[iring][ipie].fz =  zD;
	  Location[iring][ipie].fr = pow( pow(Location[iring][ipie].fx ,2.) + pow(Location[iring][ipie].fy ,2.) + pow(Location[iring][ipie].fz ,2.), 0.5);

	  float fperp =  pow( pow(Location[iring][ipie].fx ,2.) + pow(Location[iring][ipie].fy ,2.) ,0.5);

	  //cout << fperp << " " << Location[iring][ipie].fy << endl;

	  
	  Location[iring][ipie].fthetaMax = asin(fperp/Location[iring][ipie].fr);
	  Location[iring][ipie].fthetaMax = asin(fperp/Location[iring][ipie].fr);
	  Location[iring][ipie].fthetaMin = asin((fperp - dr)/Location[iring][ipie].fr);
	  Location[iring][ipie].fdtheta = Location[iring][ipie].fthetaMax -  Location[iring][ipie].fthetaMin;

	  //cout << "front thetaMax = " << Location[iring][ipie].fthetaMax << endl;
	  
	  Location[iring][ipie].fphiMin = 2.*3.1415927*(float)(ipie-Npie/2.)/(float)Npie;//atan2(Location[iring][ipie].fy, Location[iring][ipie].fx);
	  Location[iring][ipie].fphiMax = 2.*3.1415927*(float)(ipie+1 - Npie/2.)/(float)Npie;//atan2(Location[iring][ipie].by, Location[iring][ipie].bx);
	  Location[iring][ipie].fdphi = 2.*3.1415927/(float)Npie;//Location[iring][ipie].fphiMax - Location[iring][ipie].fphiMin;

	  //also need to make sure theta goes from 0 to 180
	  //also remember fragment phi goes from -180 to 180
	  //these phi are only true in the detector frame, will need to change this so we can ultimately put the detector where we want

	  
	  //cout << "iring= " << iring << endl;
	  //cout << "testing... " << Location[iring][ipie].fphiMax << endl;
	  //cout << "testing... " << fperp << endl;
	  
	}
      
    }

  //cout << "Location[nRing-1][Npie-1].fthetaMax = " << Location[nRing-1][Npie-1].fthetaMax << endl;
  
  return 1;

}

int detector::setGeometryBack(int nRing, int Npie, float rmin, float rmax)
{

  dr = (rmax-rmin)/(float)nRing;
  
  for(int ipie=0; ipie<Npie; ipie++)
    {
      for(int iring=0;iring< nRing; iring++)
	{


	  Location[iring][ipie].bx = xD;
	  Location[iring][ipie].by = yD + rmin + dr*(float)(iring+1);
	  Location[iring][ipie].bz = zD + thickness;
	  Location[iring][ipie].br = pow( pow(Location[iring][ipie].bx ,2.) + pow(Location[iring][ipie].by ,2.) + pow(Location[iring][ipie].bz ,2.), 0.5);


  	  float bperp = pow( pow(Location[iring][ipie].bx ,2.) + pow(Location[iring][ipie].by ,2.) ,0.5);

	  
	  Location[iring][ipie].bthetaMax = asin(bperp/Location[iring][ipie].br);
	  Location[iring][ipie].bthetaMin = asin((bperp - dr)/Location[iring][ipie].br);
	  Location[iring][ipie].bdtheta = Location[iring][ipie].bthetaMax -  Location[iring][ipie].bthetaMin; 


	  
	  //Location[iring][ipie].fphiMin = atan2(Location[iring][ipie].fy, Location[iring][ipie].fx);
	  //Location[iring][ipie].fphiMax = atan2(Location[iring][ipie].by, Location[iring][ipie].bx);
	  Location[iring][ipie].bphiMin = Location[iring][ipie].fphiMin;
	  Location[iring][ipie].bphiMax = Location[iring][ipie].fphiMax;
	  //Location[iring][ipie].fdphi = Location[iring][ipie].fphiMax - Location[iring][ipie].fphiMin;
	  Location[iring][ipie].bdphi = Location[iring][ipie].fdphi;

	}
    }

  //cout << "Location[nRing-1][Npie-1].bthetaMax = " << Location[nRing-1][Npie-1].bthetaMax << endl;
  
  return 1;

}



int detector::setGeometry(int nXsegments, int nYsegments, float xmin, float xmax, float ymin, float ymax, int dType)
{
  
  if(!didSetGeo)
    {
      Location = new location*[nXsegments];
      for(int ix=0; ix < nXsegments; ix++)
	{
	  Location[ix] = new location[nYsegments];
	}
      didSetGeo=1;
    }
  else
    {
      cout << endl << "Not reallocating Location matrix, if number of segments is changed beware!" << endl;
    }
  
  if(dType==0)
    {
      return(setGeometry(nXsegments,nYsegments,xmin,xmax));
      //will have to figure out better way to detector type since ymax, ymin not used
    }
  else if(dType==1) //for the rectangle (type TTT say) we'll want to compare x,y,z coordinates of fragment
    {
      return 0;
    }
  else
    {
      cout << "Functionality not added yet" << endl;
      abort();
      return 0;
    }

  return 0;
  
}


int detector::event(float theta0, float phi0, int coordinate0) //assumes particle travelling from origin (target position)
{

  float xR = ran.Rndm();
  float yR = ran.Rndm();

  // if( theta0 < Location[0][0].fthetaMin ||  theta0 > Location[Nx-1][Ny-1].fthetaMax )
  //   {
  //     return 0;
  //   }
  
  if(thickness > 0)
    {
  
      for(int ix=0; ix < Nx; ix++)
	{
	  for(int iy=0; iy < Ny; iy++)
	    {
	      if(theta0 > Location[ix][iy].fthetaMin && theta0 < Location[ix][iy].fthetaMax)
		{
		  if(theta0 > Location[ix][iy].bthetaMin && theta0 < Location[ix][iy].bthetaMax)
		    {
		      if(phi0 > Location[ix][iy].fphiMin && phi0 < Location[ix][iy].fphiMax)
			{
			  if(phi0 > Location[ix][iy].bphiMin && phi0 < Location[ix][iy].bphiMax)
			    {
			      //cout << endl << "theta0 = " << theta0 << endl;
			      thetaHit = Location[ix][iy].fthetaMin + xR*Location[ix][iy].fdtheta;
			      //cout << "thetaHit = " << thetaHit << endl;
			      phiHit = Location[ix][iy].fphiMin + yR*Location[ix][iy].fdphi;
			      segmentXhit = ix;
			      segmentYhit = iy;
			      return 1;
			    }
			}
		    }
		}
	    }
	}
    }
  else
    {
      for(int ix=0; ix < Nx; ix++)
	{
	  for(int iy=0; iy < Ny; iy++)
	    {
	      if(theta0 > Location[ix][iy].fthetaMin && theta0 < Location[ix][iy].fthetaMax)
		{
		  if(phi0 > Location[ix][iy].fphiMin && phi0 < Location[ix][iy].fphiMax)
			{
			  //cout << "theta0 = " << theta0 << endl;			  			  
			  thetaHit = Location[ix][iy].fthetaMin + xR*Location[ix][iy].fdtheta;
			  //thetaHit = atan2( (xR+(float)ix)*dr + Location[0][0].fy - dr, zD  );
			  //cout << "thetaHit = " << thetaHit << endl;
			  phiHit = Location[ix][iy].fphiMin + yR*Location[ix][iy].fdphi;
			  segmentXhit = ix;
			  segmentYhit = iy;
			  return 1;
			}
		}
	    }
	}
    }


  return 0;
  
  
}


int detector::event(float v0, float theta0, float phi0, float x0, float y0, float z0, int coordinate0)
{
  cout << "functionality will be added later" << endl;
  abort();
  return 0;

}


//*********************************************************
/**
*Destructor
*/
detector::~detector()
{
  delete Location;
}






