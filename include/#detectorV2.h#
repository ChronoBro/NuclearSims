#ifndef _detectorV2
#define _detectorV2



#include "MathGeoLib.h"
//#include "math.h"
#include <vector>
#include "TRandom.h"

/*******************************************************

What am I trying to do here?
________________________________________________________
--------------------------------------------------------

1. Create 3D detectors 
  a. Ultimately I'll want to create full 3D polyhedra
    - This could provide a problem for annular detectors

2. Ultimately want to track the DISTANCE a particle could
   travel through the detector

3. To start I'll just create flat front and back faces
   spaced by thickness and check intersection with these
  a. How will I check the total distance available in
     the detector?


*******************************************************/

class detectorV2
{

 public:
  detectorV2();
  detectorV2(float,float); //distance along z-axis and thickness
  detectorV2(float,float,float,float); //x,y,z position and thickness
  detectorV2( const detectorV2 &obj); //copy constructor
  ~detectorV2();
  
  int setGeometry(int type0, int nxstrips0, int nystrips0, float x1min0, float x1max0, float x2min0, float x2max0);
                             //preload geometry, I'll start with front and back rectangular polygons, as well as annular devices
                             //type=1 rectangular type=2 circular, if type=0 geometry is not set
                             //need to set min/max of plane dimensions, will handle these differently depending on type
                             //this will assume back face is the same as front face at first
  
  int setGeometryBack(int type0, int nxstrips0, int nystrips0, float x1min0, float x1max0, float x2min0, float x2max0);
                             //explicitely change the back face geometry
  int setDetNormal();
  int setDetFrame(); //used to create LOCAL coordinates on the detector, strips will be created in this frame
  int event(Ray trajectory);
  int setNormal(vec);
  int setLocalBasis();
  void assignMeasurement(float x, float y); //using 2D LOCAL coordinates 
  void assignMeasurement(Ray trajectory); //just for solid AABB (rectangles aligned with world-axes)
  
  float distTravelled;
  float distInteraction;
  float measuredDist;
  Ray measuredTrajectory;

  
 private:

  Circle annularFront;
  Circle annularBack;
  Plane detFrontPlane;
  Plane detBackPlane;
  AABB aabbDet;
  OBB obbDet;
  
  float3 pointIntersect;
  TRandom ran;
 
  
  float thickness;
  
  vec center;
  vec normal;
  
  vec detXaxis;
  vec detYaxis;

  int type;
  int nxstrips;
  int nystrips;
  float x1min;
  float x1max;
  float x2min;
  float x2max;
  

};

#endif
