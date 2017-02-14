#include <string>
#include <fstream>
#include <iostream>

using namespace std;

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

  CLoss(string,float);
  ~CLoss();

    float getEout(float,float);
    float getEin(float,float);
    float getDedx(float);
};
