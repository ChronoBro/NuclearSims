#include "plane.h"

using namespace std;

plane::plane()
{
  //printf("%s\n", "you may now use plane functions but remember to set plane normal and center");
  didSet = 0;
}


plane::plane(double * n0, double * p0) //p will be plane center, normX & normY will be assummed to be from this point
{
  didSet = 1;
  for (int i=0;i<3;i++)
    {
      n[i] = n0[i];
      p[i] = p0[i];
    }

}

void plane::findNormal(double *p10, double *p20, double *p30)
{

  double cross1[3];
  double cross2[3];

  //finding normal vector via cross product
  for (int i=0;i<3;i++)
    {
      cross1[i] = p20[i] - p10[i];
      cross2[i] = p30[i] - p10[i];
    }
  
  n[0] = cross1[1]*cross2[2] - cross1[2]*cross2[1];
  n[1] = cross1[2]*cross2[0] - cross1[0]*cross2[2];
  n[2] = cross1[0]*cross2[1] - cross1[1]*cross2[0];

  double mag=0.;

  for (int i=0;i<3;i++)
    {
      mag += pow(n[i],2.);
    }

  mag = sqrt(mag);

  //creating unit vector
  for(int i=0;i<3;i++)
    {
      n[i] = n[i]/mag;
    }

  //check that it's pointing towards center
  double dot;

  dot = n[0]*p[0] + n[1]*p[1] + n[2]*p[2];
  if(dot > 0)
    {
      for(int i=0;i<3;i++)
	{
	  n[i] = -1.*n[i];
	}
    }

}

void plane::setPlaneCenter(double *p0)
{
  for(int i=0;i<3;i++)
    {
      p[i] = p0[i];
    }
}

void plane::setNormal(double *n0)
{
  //check if actually a normal vector...
  double mag=0.;
  
  for(int i=0;i<3;i++)
    {
      mag += pow(n[i],2);
      n[i] = n0[i];
    }

  mag = sqrt(mag);
  if(mag!=1)
    {
      cout << "did not input a normalized vector! It will be rescaled." << endl;

      for(int i=0;i<3;i++)
	{
	  n[i] = n0[i]/mag;
	}
      
    }
  
  
}

//only works for cartesian vectors
double * plane::planeCoordinate(double * p0)
{
  planeVector[0] = normX[0]*p0[0] + normX[1]*p0[1] + normX[2]*p0[2];
  planeVector[1] = normY[0]*p0[0] + normY[1]*p0[1] + normY[2]*p0[2];

  return planeVector;
}


double * plane::labCoordinate(double xprime, double yprime)
{

  for(int i=0;i<3;i++)
    {
      labVector[i] = xprime*normX[i] + yprime*normY[i] + p[i];
    }

  return labVector;

}



int plane::setPlaneCoordinates(double * normX0, double * normY0)
{

  for (int i = 0;i<3;i++)
    {
      normX[i] = normX0[i] - p[i];
      normY[i] = normY0[i] - p[i];
    }

  return 1;
  
}



void plane::checkSet()
{
  if(!didSet)
    {
      cout << "you have to set normal vector and center point! Aborting... " << endl;
    }
  abort();
  
}

void plane::moveCenter(double x0, double y0, double z0)
{
  p[0]=x0;
  p[1]=y0;
  p[2]=z0;
}


//returns 0 if parallel, -1 if IN PLANE, otherwise returns paramter t that gives intersection
double plane::lineIntersect(double * p0, double * v0) //position vector and velocity vector of line
{
  double dot1=0.;
  double dot2=0.;
  double diff[3];

  for(int i=0;i<3;i++)
    {
      diff[i] = p[i]-p0[i];
    }

  dot1 = n[0]*diff[0]+n[1]*diff[1] + n[2]*diff[2];

  if(dot1 == 0)
    {
      return 0;
    }

  dot2 = n[0]*v0[0]+n[1]*v0[1]+n[2]*v0[2];

  if(dot2 > 0.000001)
    {
      return dot1/dot2;
    }
  else
    {
      return -1.;
    }

}


plane::~plane()
{


}
