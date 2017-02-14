#include "detector.h"


/**
 * Constructor
 \param Z0 is the atomic number of fragment
 \param mass0 is mass of fragment in amu
 \param filename specifies the name of the file with energy loss information
 \param CsI_res0 fractional energy resolution of this fragment in CsI
 \param thickness is target thickness in mg/cm2 
*/


detector::detector(float dist0, float thickness0)
{
  dist = dist0;
  thickness = thickness0;
  
}

detector::detector(float x0, float y0, float z0, float thickness0)
{
  x=x0;
  y=y0;
  z=z0;
  thickness = thickness0;

}

int detector::setGeometry(int nRing, int Npie, float rmin, float rmax)
{

  xy = new float*[nRing];

  for(int iring=0; iring < nRing; iring++)
    {
      xy[iring] = new float[Npie];
    }

  float dr = (rmax-rmin)/(float)nRing;
  
  for(int ipie=0; ipie<Npie; ipie++)
    {
      for(int iring=0;iring< nRing; iring++)
	{
	  
	}
    }

  return 0;
}

int detector::setGeometry(int nXsegments, int nYsegments, float xmin, float xmax, float ymin, float ymax, int dType)
{

  if(dType==0)
    {
      return(setGeometry(nXsegments,nYsegments,xmin,xmax));
      //will have to figure out better way to detector type since ymax, ymin not used
    }
  else if(dType==1)
    {
      return 0;
    }
  else
    {
      cout << "Functionality not added yet" << endl;
      return 0;
    }

  return 0;
  
}


//*********************************************************
/**
*Destructor
*/
detector::~detector()
{
  delete xy;
}






