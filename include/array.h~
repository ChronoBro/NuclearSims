
#ifndef _array
#define _array

#include "tower.h"
#include <string>
#include <fstream>

using namespace std;

/**
 *!\brief handle the Hira array
 *
 * contains the geometrical properties of the HiRA array and 
 * indicates is a fragment is detected in the Si detector
 */

class CArray
{
 public:
  CTower * Tower[5]; //!< parameter of each tower 
  CArray(float dist0, float fromCenter0, float Active, float Physical);
  ~CArray();
  int hit(float,float,float,float);
  int hitTower; //!< number of tower in which fragment was detected
  float pi;   //!< 3.14159...
  int nTower; //!< number of towers in the array


  float thetaRecon; //!< reconstructed theta angle in radians of detected fragment
  float phiRecon;   //!< reconstructed phi angle in radians

  float alphaCenter1; 
  float deltaAlphaTower;
  float deltaBetaTel;
  float fromCenter;
  float extra;
  float xActive;
  float yActive;
  float dist;



};


#endif
