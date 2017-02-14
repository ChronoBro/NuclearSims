#include "ttt.h"

CRandom ttt::ran;

ttt::ttt(float dist0,float theta_D0, float height0, float width0,float offset0, int N_v0, int N_h0)
{
  offset = offset0;
  dist = dist0;
  theta_D = theta_D0;
  height = height0;
  width = width0;
  N_v = N_v0;
  N_h = N_h0;
  chi = dist*(cos(theta_D)+pow(sin(theta_D),2)/cos(theta_D));
  dl = sqrt(1+pow(tan(theta_D),2)); // not actual length of line/detector_plane, just transform between distance along line to x distance so this quantity is unitless
  dy = height/N_h;
  dx = (width/N_v)/dl;
  y_max = height/2+offset;
  y_min = -height/2+offset;
  x_2 = dist*sin(theta_D) + (width/2)/dl; //x of center plane = dist*sin(theta_D) so these equations tell the x extrema of the detector
  x_1 = dist*sin(theta_D) - (width/2)/dl;
}

int ttt::hit(float theta, float phi)
{
  
  float yr = ran.Rndm();
  float xr = ran.Rndm();
  hitx=0;
  hity=0;
  if (x_1 < x_2)
    {x_min = x_1;
     x_max = x_2;
    }
  else
    {x_min = x_2;
     x_max = x_1;
    }

  //Calculting (theta,phi) projection on detector plane

  a = pow(sin(theta),2)*(1+pow(cos(phi),2)*pow(tan(theta_D),2))-1;
  b = -2*sin(theta)*cos(phi)*tan(theta_D)*chi;
  c = pow(chi,2);

  //these values are what you get solving for rho with rho^2=x^2+y^2+z^2 and the defined plane equation

  //  if ((b*b-4*a*c) < 0) 
  //  return (0);

  
  rho1 = (-b + sqrt(b*b-4*a*c))/(2*a);
  rho2 = (-b - sqrt(b*b-4*a*c))/(2*a);
  
  // cout << "rho1=" << rho1 << endl;
  // cout << "rho2=" << rho2 << endl;


  if (rho1 <dist)
    rho = rho2;
  else 
    rho = rho1;

  //  cout << "rho=" << rho << endl;
  //this works as the sphere rho=dist is tangent to the plane so rho along the plane will always be greater than rho=dist //might need to fix if both solutions are greater at large angles

  //(x,y,f(x)) mapping to (rho(theta,phi),theta,phi) along detector plane
  x = rho*sin(theta)*cos(phi);
  y = rho*sin(theta)*sin(phi); 
  z = rho*cos(theta);
  
 
 
  // if((x*x+y*y) > (pow(tan(acos(-1)/25),2)*z*z)) return (0); //this excludes the longitudinal/transverse (</>) events with a cone of fixed theta 

  //cout << "x="<<x<<", x_max="<<x_max<<", x_min=" << x_min<<endl;
   // cout << "y="<<y<<", y_max="<<y_max<<", y_min="<<y_min<<"\n"<<endl;
  //logical telling if detector was hit and specifying section hit
  if((x_min<x && x<x_max)&&(y_min<y && y<y_max))
    {
      hitx = (int)(x-x_min)/dx;
      hity = (int) (y-y_min)/dy;
      // cout << "hitx= " << hitx << endl;
      // cout << "hity= " << hity << endl;
      x = dx*(xr+(float)hitx)+x_min;
      y = dy*(yr+(float)hity)+y_min;
      thetaHit = atan2(sqrt(x*x+y*y),z);
      phiHit = atan2(y,x);
      CsIstrip1 = (int) (x-x_min)/dx*6/N_v;
      CsIstrip2 = (int) (y-y_min)/dy*6/N_h;
     // cout << "CsIx= " << CsIstrip1 << endl;
     // cout << "CsIy= " << CsIstrip2 << endl;
     // cout << "theta="<<theta<<","<<thetaHit << endl;
      // cout << "phi="<<phi<<","<<phiHit<<"\n"<<endl;
      return 1;

    }
  
  return 0;
}
