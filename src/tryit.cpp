#include "loss.h"
#include <iostream>
#include <string>

using namespace std;

int main()
{
  string filename("deuteron.loss");
    loss Loss(filename.c_str());

  float energy = 41.5197;
  float thick = 76.3;
  cout << Loss.GetEout(energy,thick) << endl;

}
