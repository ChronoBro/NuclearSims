#ifndef _tele
#define _tele



#include <iostream>
#include "random.h"
#include "sle.h"

/**
 *!\brief  deals with the parameters of a single HiRA telescope
 *
 * stores location and indicates if a particle hit a HiRA telescope
 */
class CTele
{
public:
  static CRandom ran;
  float pi;
  float r_center[3]; //!<coordinates of center of the front face (cm)
  //float r_normal[3]; //!<unit vector perpendicular to the front face (cm)
  float r_front[3]; //!<unit vector pointing along front strips (cm)
  float r_back[3]; //!<unit vector pointing along back strips (cm)

  float xActive; //!< in-plane dimension of active region of Si detector (cm)
  float yActive; //!< out-of-plane dimension

  CTele(float*,float*,float*,float,float);
  CTele(){};
  void init(float*,float*,float*,float,float);


  int hit(float,float,float,float); 

  int ixStrip; //!< in-plane strip number that particle hit
  int iyStrip; //!< out-of-plane strip number
  int ICsI; //!< id of CsI crystal which was hit
  int ICsI3;

  float thetaRecon,phiRecon;



};

#endif
