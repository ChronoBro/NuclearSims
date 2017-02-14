#include "polyScat.h"

polyScat::polyScat(float thickness)
{
  //thickness is in mg/cm2

  Zhydrogen = 1.;
  Zcarbon = 6.;
  Nhydrogen =  thickness/1000./14.*2.*6.02e23; 
  Ncarbon = thickness/1000./14.*6.02e23;
  Hydrogen = new multScat(Nhydrogen,Zhydrogen);
  Carbon = new multScat(Ncarbon,Zcarbon);
}
//***********************************************
float polyScat::thetaRMS(float Zproj, float Eproj,float fractionalThick)
{
  float HydrogenRMS = Hydrogen->thetaRMS(Zproj,Eproj,fractionalThick);
  float CarbonRMS = Carbon->thetaRMS(Zproj,Eproj,fractionalThick);
  float total = sqrt(pow(HydrogenRMS,2)+pow(CarbonRMS,2));
  return total;
}
