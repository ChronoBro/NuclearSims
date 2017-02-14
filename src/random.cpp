#include "random.h"

TRandom3 CRandom::random(33456);


double CRandom::Rndm()
{
  return  random.Rndm();
}
//***********************
float CRandom::Gaus(float mean, float sigma)
{
  // returns random number with a gaussian distribution
  //with given mean and stabdard deviation (sigma)
  return random.Gaus(mean,sigma);
}
//*************************
float CRandom::expDecayTime(float width)
{
  // returns a time from an exponential decay distribution
  //consistent with the total decay width
return - 0.65824*log(Rndm()+1.0e-37)/width;
}
//****************************************
float CRandom::BreitWigner(float mean, float width)
{
  return random.BreitWigner(mean,width);
}
