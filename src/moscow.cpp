#include "moscow.h"


const float moscow::constant= sqrt(979.82);

moscow::moscow(string name)
{
  ifFile.open(name.c_str());
  if (ifFile.fail())
    {
      cout << "file not found " << name << endl;
    }

  ifFile >> ntot >> A1 >> A2 >> A3 ;
}
bool moscow::getEvent(float Etotal,CFrame * Frag1,CFrame *Frag2,CFrame *Frag3)
{
  float Etot;
  float p1[3];
  float p2[3];
  float p3[3];

     ifFile >> Etot >> p1[0] >> p1[1] >> p1[2] >> p2[0] >> p2[1] >> p2[2] >>
    p3[0] >> p3[1] >> p3[2];

     if (ifFile.eof()) abort();
     if (ifFile.bad()) abort();


    float tot=0;
    for (int i=0;i<3;i++) 
      {
        Frag1->v[i] = p1[i]/A1/constant*sqrt(Etotal/Etot);
        Frag2->v[i] = p2[i]/A2/constant*sqrt(Etotal/Etot);
        Frag3->v[i] = p3[i]/A3/constant*sqrt(Etotal/Etot);
      }

    return 1;
}

bool moscow::getEvent8C(float Etotal,CFrame * Frag1,CFrame *Frag2,CFrame *Frag3)
{
  float Etot;
  float p1[3];
  float p2[3];
  float p3[3];

     ifFile >> Etot >> p1[0] >> p1[1] >> p1[2] >> p2[0] >> p2[1] >> p2[2] >>
    p3[0] >> p3[1] >> p3[2];

     if (ifFile.eof()) abort();
     if (ifFile.bad()) abort();


    float tot=0;
    for (int i=0;i<3;i++) 
      {
        Frag1->v[i] = p1[i]/Frag1->A/constant;
        Frag2->v[i] = p2[i]/Frag2->A/constant;
        Frag3->v[i] = p3[i]/Frag3->A/constant;
        tot += pow(Frag1->v[i]/.9784,2)/2.*Frag1->A;
        tot += pow(Frag2->v[i]/.9784,2)/2.*Frag2->A;
        tot += pow(Frag3->v[i]/.9784,2)/2.*Frag3->A;
      }

    for (int i=0;i<3;i++) 
      {
        Frag1->v[i] *= sqrt(tot/Etotal);
        Frag2->v[i] = p2[i]/(tot/Etotal);
        Frag3->v[i] = p3[i]/(tot/Etotal);
      }

    return 1;
}

