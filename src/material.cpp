#include "material.h"

using namespace std;

material::material()
{

}

material::material(string name0, double density0, double Z0, double A0)
{
  setMaterial(name0, density0, Z0, A0);
}





int material::setMaterial(string &name0, double &density0, double &Z0, double &A0)
{
  density = density0;
  Z = Z0;
  A = A0;
  name = name0;

  string filename = "../loss/";
  filename.append(name.c_str());
  
  string copy = filename;
  copy.append("/Hydrogen.loss");
  lossH.setParameters(copy); // from here out I have to add a mass argument to the ELost

  copy = filename;
  copy.append("/Helium.loss");
  lossHe.setParameters(copy);

  copy = filename;
  copy.append("/Lithiumloss");
  lossLi.setParameters(copy);

  copy = filename;
  copy.append("/Calcium.loss");
  lossCa.setParameters(copy);

  copy = filename;
  copy.append("/Argon.loss");
  lossAr.setParameters(copy);
  
  
  return 1;
}

double material::Elost(float energy0, double length, int particleZ0, double particleA0)
{
  

  if(particleZ0==1)
    {
      return lossH.getEout(energy0 , length , particleA0);
    }
  else if(particleZ0==2)
    {
      return lossHe.getEout(energy0, length, particleA0);
    }
  else if(particleZ0==3)
    {
      return lossLi.getEout(energy0, length, particleA0);
    }
  else if(particleZ0==18)
    {
      return lossAr.getEout(energy0, length, particleA0);
    }
  else if(particleZ0==20)
    {
      return lossCa.getEout(energy0, length, particleA0);
    }
  else
    {
      cout << "Element Z = " << particleZ0 << " not currently supported. Go ahead and add the loss files (for each material!) and add it to material.cpp!" << endl;
      return 0.;
    }
  
}



material::~material()
{

}
