#include "planeSegment.h"

//planes are defined by normal vector and then a point on the plane
//to get a proper plane segment one needs the boundaries of the plane 


planeSegment::planeSegment(const plane &plane0, int nBP)
{

  pl = &plane0;
  NBP = nBP;

  for(int i=0;i<NBP;i++)
    {
      boundaryPoints[i] = new double[3];

    }
  
}




planeSegment::~planeSegment()
{
  for(int i=0;i<NBP;i++)
    {
      delete boundaryPoints[i];
    }
  

}
