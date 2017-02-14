#include "tower.h"

//**************************************************************
  /**
   * Constructor to set up a tower of HiRA detectors
   \param nTele0 is number of telescopes in tower
   \param rad_si is radius from nominal target location in cm to the front of the Si detectors
   \param delta_dist is the real target location relative to the nomimal nominal in cm. Positive values are further away from the Hira array.
   \param alpha is the angle in the reaction plane of the center of the tower relative to the beam axis
   \param beta0 is an array  of out-of-plane angles for each telescope.
   \param xActive0 is the active dimension of the Si telescope in plane (cm)
   \param yActive0 is the active dimension of the Si telescope out-of-plane (cm)
   */
CTower::CTower(int nTele0,float rad_si,float delta_dist,float alpha, 
float* beta0, float xActive0, float yActive0)
{
  nTele = nTele0; 
  betaA = new float [nTele];
  for(int i=0;i<nTele;i++) betaA[i] = beta0[i];
  Tele = new CTele * [nTele];
  pi = acos(-1.);
  xActive = xActive0;
  yActive = yActive0;
  for (int itele =0;itele<nTele;itele++)
    {
      //out of plane angle
      float beta = betaA[itele];

      //theta and phi angle of telescope
      float theta = acos(cos(beta)*cos(alpha));
      float phi = asin(beta/theta);
      if (alpha < 0.) phi = pi - phi;
      
      // xyx coordinates of center of telescopes 
      float r_center[3];
      r_center[0] = rad_si*sin(theta)*cos(phi);
      r_center[1] = rad_si*sin(theta)*sin(phi);
      r_center[2] = rad_si*cos(theta) + delta_dist;

      /*
      // xyx coordinates of normal to center
      float r_normal[3];
      r_normal[0] = -sin(theta)*cos(phi);
      r_normal[1] = -sin(theta)*sin(phi);
      r_normal[2] = -cos(theta);
      */

      //xyz coordinates of normal vector perpt to front strips
      float r_front[3];
      r_front[0] = cos(alpha);
      r_front[1] = 0.;
      r_front[2] = -sin(alpha);


      //xyz coordinates of normal vector perpt to back strips
      float r_back[3];
      r_back[0] = sin(alpha)*sin(beta);
      r_back[1] = -cos(beta);
      r_back[2] = cos(alpha)*sin(beta);


      Tele[itele] = new CTele (r_center,r_front,r_back,xActive,yActive);
    }
}
//*****************************************************
  /**
   * Destructor 
   */
CTower::~CTower()
{
  for (int i=0;i<nTele;i++) delete Tele[i];
  delete [] Tele;
  delete [] betaA; 
}
//***********************************************************
  /**
   * yes yes or no to question of whether a particle hits a Si detector 
   * in the tower
   \param theta is the zenith angle of particle in radians.
   \param phi is the polar angle of particle in radians.
   \param xtarget is the in-plane offset of interaction in target from nomimal beam center (cm)
\param ytarget is the out-of-plane offset
   */
int CTower::hit(float theta,float phi, float xtarget, float ytarget)
{
  int nhit = 0;
  for (int itele=0;itele<nTele;itele++)
    {
      nhit += Tele[itele]->hit(theta,phi,xtarget,ytarget);
      if (nhit == 1) 
	{
	  hitTele = itele;
          thetaRecon = Tele[itele]->thetaRecon;
          phiRecon = Tele[itele]->phiRecon;
          break;
	}
    }
  return nhit;
}
