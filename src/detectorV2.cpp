#include "detectorV2.h"

//MATH_BEGIN_NAMESPACE

using namespace math;

detectorV2::detectorV2(const detectorV2 &obj)
{
  type = obj.type;
  center = obj.center;
  thickness = obj.thickness;
  normal = obj.normal;
  setLocalBasis();
  detFrontPlane.normal = normal;
  detBackPlane.normal = normal;
  
  detFrontPlane.d = center.Length();
  detBackPlane.d = center.Length() + thickness;
}

detectorV2::detectorV2()
{
  type=0;
}

detectorV2::detectorV2(float z0, float thick)
{
  type=0;
  center.Set(0.,0.,z0);
  thickness = thick;

  //by default I'll have it pointing towards target, can set normal later
  normal = center.Normalized();
  for(int i=0;i<3;i++)
    {
      normal[i] = -1.*normal[i];
    }
  setLocalBasis();

  //since only the plane implementation of this library will give me the point of intersection I can check if trajectory
  //intersects the object and then have the plane return the point of intersection
  
  detFrontPlane.normal = normal;
  detBackPlane.normal = normal;

  detFrontPlane.d = center.Length();
  detBackPlane.d = center.Length()+thickness;
}


detectorV2::detectorV2(float x0, float y0, float z0, float thick)
{
  type=0;
  center = float3(x0,y0,z0);
  thickness = thick;

  //by default I'll have it pointing towards target, can set normal later
  normal = center.Normalized();
  for(int i=0;i<3;i++)
    {
      normal[i] = -1.*normal[i];
    }
  setLocalBasis();

  detFrontPlane.normal = normal;
  detBackPlane.normal = normal;

  detFrontPlane.d = center.Length();
  detBackPlane.d = center.Length()+thickness;
  //std::cout << "detectorV2(x,y,z,thickness) constructor called" << std::endl;
}



detectorV2::~detectorV2()
{


}

int detectorV2::setNormal(vec direction)
{
  normal = direction;
  if(normal.Length()!=1)
    {
      normal.Normalize();
    }
  setLocalBasis();

  
  return 1;
  
}


int detectorV2::setLocalBasis()
{
  detXaxis = normal.Perpendicular();
  detYaxis = normal.AnotherPerpendicular();
  
  return 1;
}



int detectorV2::setGeometry(int type0, int nxstrips0, int nystrips0, float x1min0, float x1max0, float x2min0, float x2max0)
{
  type=type0;
  nxstrips=nxstrips0;
  nystrips=nystrips0;
  x1min = x1min0;
  x2min = x2min0;
  x1max = x1max0;
  x2max = x2max0;

  
  if(!type)
    {
      std::cout << "Detector Geometry not set!" << std::endl;
      return 0;
    } 
  else if(type==1)
    {
      //rectangular shape
      //Remember rectangular dets currently only oriented along world axes
      float3 size;
      size[0]=x1max-x1min;
      size[1]=x2max-x2min;
      size[2]=thickness;

      
      float3 centerPoint = center + center.unitZ.ScaledToLength(thickness/2.);

      aabbDet.SetFromCenterAndSize(centerPoint,size);
      
    }
  else if(type==2)
    {
      //annular shape
      annularFront.pos = center;
      float3 offset = normal;

      for(int i=0;i<3;i++)
	{
	  offset[i] = -1*offset[i]*thickness;
	}
      
      annularBack.pos = center + offset;
      annularFront.normal = normal;
      annularFront.normal = normal;
      
    }

  return 1;
  
}


int detectorV2::event(Ray trajectory)
{
  float4 d = trajectory.dir.ToPos4();
  float * distIntersect;
  //*distIntersect = 0.;
  
  
  if(type==0)
    {
      std::cout << "Geometry not set to take events! Aborting..." << std::endl;
      abort();
    }
  else if(type==1)
    {
      //std::cout << "Functionality not implemented yet" << std::endl;
      float dnear;
      float dfar;
      if(trajectory.Intersects(aabbDet, dnear, dfar))
	{
	  distTravelled = dfar - dnear;
	  distInteraction = dnear + ran.Rndm()*distTravelled;
	  assignMeasurement(trajectory);
	  return 1;
	}
      else
	{
	  return 0;
	}
    }
  else if(type==2)
    {
      if(trajectory.IntersectsDisc(annularFront))
	{	  
	  trajectory.Intersects(detFrontPlane, distIntersect);
	  pointIntersect = trajectory.GetPoint(*distIntersect);
	  if(trajectory.IntersectsDisc(annularBack))
	    {
	      //if(event) then make function to calculate 'measured' trajectory
	      return 1;
	    }
	  else
	    {
	      
	    }
	}
      else
	{
	  return 0;
	}
      
    }
  
  return 0;
  
}

void detectorV2::assignMeasurement(Ray trajectory)
{
  float zRes = 2.; //assuming z-resolution FWHM in inches
  measuredTrajectory.dir = aabbDet.PointInside(ran.Rndm(), ran.Rndm(), ran.Rndm()).Normalized();
  float3 realD = trajectory.GetPoint(distInteraction);
  float3 reconD = realD;
  
  for(;;)
    {
      reconD[2] = realD[2] + ran.Gaus(0,zRes*2.54/2.355);
      measuredDist = reconD.Length();
      measuredTrajectory.dir[2] = reconD[2];
      if(measuredTrajectory.dir[2] > aabbDet.MinZ() || measuredTrajectory.dir[2] < (aabbDet.MinZ()+ thickness) )
	{
	  break;
	}
    }
  measuredTrajectory.dir.Normalize();
  float3 zero = zero.zero;
  measuredTrajectory.pos = zero;
}


void detectorV2::assignMeasurement(float x, float y)
{
  if(type==1)
    {
      //will make this rectangular shit easy for now 
    }
  else if(type==2)
    {


    } 
  
}
