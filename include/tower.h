#ifndef _tower
#define _tower

#include "tele.h"



/**
 *!\brief handles the setup of a tower in the HiRA array
 *
 * This structure defines the parameters of a tower and allows one
 * to determine if a particle is detected in the tower
 */

class CTower
{
 public:
 CTele **Tele; 
 CTower(int,float,float,float,float*,float,float);
 ~CTower();
 int nTele; //!< number of telescopes in tower
 float * betaA; //!< beta angles (out-of-plane angle) of the telescopes in the tower
 float alpha;  //!< in-plane  angle of tower [rad]
 float xActive; //!< in-plane actice dimension of Si detector
 float yActive; //!< out-of-plane active dimension of Si detector 

 int hit(float,float,float,float);
 int hitTele;//!< specifies which telescope number  was hit
 float pi; //!< 3.14159....

 float thetaRecon; //!< specifies the reconstructed theta angle of detected fragment
 float phiRecon; //!< specifies the reconstructed phi angle [rad]
};


#endif
