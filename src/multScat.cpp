#include "multScat.h"

multScat::multScat(float totalThick0,float Ztarget0):pi(acos(-1.)),
				     alpha(1./137.), e2(1.44e-13)
{
  Ztarget = Ztarget0;
  totalThick = totalThick0;
  Lrad = log(184.15/pow(Ztarget,(float)(1./3.)));
  LdashRad = log(1194./pow(Ztarget,(float)(2./3.)));
  f = pow(Ztarget*alpha,2);
  f = f*(1./(1+f) + 0.20206 - 0.0369*f + .0083*pow(f,2));

  factor = sqrt(4*pi*totalThick*(pow(Ztarget,2)*(Lrad-f)+Ztarget*LdashRad))
            *e2; 
}

//*********************************************
   float multScat::thetaRMS(float Zproj, float Eproj,float fractionalThick)
{
  return factor*sqrt(fractionalThick)*Zproj/Eproj;
}
