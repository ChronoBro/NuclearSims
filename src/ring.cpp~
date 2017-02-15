#include "ring.h"


CRandom ring::ran;  //random number generator

//constructor
ring::ring(float dist0, float rmin0, float rmax0, int Nring0, int Npie0, int place0)
{
  dist = dist0;  // distance of ring counter to target
  rmin = rmin0;  // minimum active radius of ring counter
  rmax = rmax0;  // maximum active radius of ring counter
  Nring = Nring0; // number of annular ring sectors
  Npie  = Npie0;  //number of pie sectors
  place = place0;
  theta_min = atan2(rmin,dist);
  theta_max = atan2(rmax,dist);
  dr = (rmax-rmin)/(float)Nring;
  dtheta = (theta_max - theta_min)/(float)Nring;
  dphi = 2*acos(-1)/(float)Npie;

  for(int i = 0; i < Nring; i++){
    theta_min_vec.push_back(atan2(i*dr + rmin,dist));
    theta_max_vec.push_back(atan2((i+1)*dr + rmin,dist));
  }
}

// returns 1 if particle hits detector, else returns zero
// inputs  theta = emission zenith angle in radians of scattered particle
//         phi =   emission azimuth angle in radians of scattered particle
//         phi goes from -pi to pi
int ring::hit(float theta, float phi)
{
  //  cout << "now phi is " << phi << endl;
  float x = ran.Rndm();
  float y = ran.Rndm();
  // float dead = 1.0/(2*acos(-1)*rmax);
  // if (y<dead) return (0); // getting rid of section between pies, just input "angular difference"/360 //integer division gives int back
  // if (y>(1-dead)) return (0);
  //specify ring and pie section which were hit
  if (theta_min<theta && theta<theta_max)
    {
	
      for (int j=0; j<Nring; j++)
	{
	  if(theta_max_vec.at(j) > theta && theta > theta_min_vec.at(j))
	    { hitRing=j;break;}
	}	    
      if(phi>=0)
	hitPie = int(phi/dphi);
      else
	hitPie = int(phi/dphi)-1;
      thetaHit = atan2((x+(float)hitRing)*dr+rmin,dist);
      //thetaHit = dtheta*(x+(float)hitRing)+theta_min;
      phiHit = dphi*(y+(float)hitPie);
      if (place == 0)
	{
	  segmentation_ring = 64;  //front
	}

      else if (place == 1)
	{
	  segmentation_ring = 65; //back
	}

      if(hitRing >= segmentation_ring){
	CsIstrip = (int) hitPie*16/Npie;
	CsIring = (int) 0;
      }
      else{
	CsIstrip = (int) hitPie*4/Npie;  //inner ring
	CsIring = (int) 1;
      }
      // if (phi<0. || phi > 2.*acos(-1))
// 	{
// 	  cout << endl;
// 	  cout << "phi = " << phi << endl;
// 	  cout << "theta = " << theta << endl; //for testing
// 	  cout << "thetaHit = " << thetaHit << endl;
// 	  cout << "# of Rings = " << Nring << endl;
// 	  cout << "hitRing = " << hitRing << endl;
// 	  cout << "hitCsIstrip = " << CsIstrip << endl;
// 	  cout << "hitCsIring = " << CsIring << endl;
// 	}
//      cout << "finally phi is " << phi << endl;
      please=1;
      return 1;
    }
  else
    please = 0;
    return (0);
}
