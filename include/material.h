#ifndef _material
#define _material

#include <string>
#include "math.h"
#include <iostream>
#include "loss.h"



class material
{
 public:
  
  material();
  material(std::string name0);
  material(std::string name0, double density0, double Z0, double A0);
  ~material();
  
  int setMaterial(std::string &name0, double &density0, double &Z0, double &A0);
  int setMaterial(std::string name0);
  double Elost(float energy0, double length, int particleZ0, double particleA0);
  
 private:

  double density;
  std::string name;
  double Z;
  double A;
  CLoss lossH;
  CLoss lossHe;
  CLoss lossLi;
  CLoss lossBe;
  CLoss lossCa;
  CLoss lossAr;
  
};
  
#endif
