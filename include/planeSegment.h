#ifndef _planeSegment
#define _planeSegment

#include "math.h"
#include "plane.h"






class planeSegment
{
 public:

  planeSegment(const plane& plane0, int nBP);
  ~planeSegment();
  //void setBoundary(double * p10, double * p20, double * p30, double * p40);
  
 private:

  double ** boundaryPoints;
  const plane * pl;
  int NBP;
};



#endif
