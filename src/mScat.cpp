#include "mScat.h"


float const CMultScat::a0 = 0.529e-8;
float const CMultScat::pi = 3.14159;

/**
 * Constructor
 \param Zpro0 is particle's atomic number
 \param Ztar0 is target's atomic number
 \param totalThick0 is targts thickness in mg/cm2
 */

CMultScat::CMultScat(int Zpro0, int Ztar0, float totalThick0)
{
  Zpro = Zpro0;
  Ztar = Ztar0;
  totalThick = totalThick0;  //atoms per centermeter sqaured

  a = 0.885*a0/sqrt(pow((double)Zpro,2./3.)+pow((double)Ztar,2./3.));

  totalTau = pi*pow(a,2)*totalThick;

  factor = 16.26/(float)(Zpro*Ztar)/
    sqrt(pow((double)Zpro,2./3.)+pow((double)Ztar,2./3.))*1000.;

  //cout << a << " " << totalTau << endl;
}

/**
 * return rms scattering angle 
  \param energy is the particles energy in MeV
  \param is the fractional thickness of the target through which the particle traverses 
 */

float CMultScat::thetaRMS(float energy, float fractionalThickness)
{
  float tau = fractionalThickness*totalTau;
  float alphaBar = pow((double)tau,.55);
  //cout << "alphaBAr= " << alphaBar << " factor " << factor << endl;
  float alpha = alphaBar/energy/factor; //half width at half maximum in radians
  return alpha*.851;
}
