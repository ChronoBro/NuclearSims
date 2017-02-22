#include <string>
#include <fstream>
#include <iostream>


/**
 * !\brief energy loss of particles in an absorber
 */

class CLoss
{
 public:
  int N;
  float *Ein;
  float *dedx;
  float mass;

  
  CLoss();
  CLoss(std::string,float);
  ~CLoss();
  int setParameters(std::string &);
  
    float getEout(float,float);
    float getEout(float,float,float);
    float getEin(float,float);
    float getDedx(float);
};
