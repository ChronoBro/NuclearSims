


#include "CMomDist.h"


CRandom CMomDist::random;
int const CMomDist::n=221;

/**
 * Constructor
 */
CMomDist::CMomDist()
{
  //cout << endl << "creating momentem distributions from MomDist" << endl;

  xArrayz = new float [n];
  yArrayz = new float [n];
  xArraytr = new float [n];
  yArraytr = new float [n];
  ifstream ifFilez("../../MomDist/zn02_allM.out");
  if (ifFilez.fail()) 
    {
      cout << " could not open zn02_allM.out " << endl;
      abort();
    }

  for (int i=0;i<n;i++)
    {
      ifFilez >> xArrayz[i] >> yArrayz[i];
      if (ifFilez.eof()) abort();
      if (ifFilez.bad()) abort();
      if ( i > 0) yArrayz[i] += yArrayz[i-1];
    }

  for (int i=0;i<n;i++) yArrayz[i] /= yArrayz[n-1];
  ifFilez.close();

  ifstream ifFiletr("../../MomDist/trn02_allM.out");
  if (ifFiletr.fail())
    {
      cout << " could not open trn02_allM.out " << endl;
      abort();
    }

  for (int j=0;j<n;j++)
    {
      ifFiletr >> xArraytr[j] >> yArraytr[j];
      if (ifFiletr.eof()) abort();
      if (ifFiletr.bad()) abort();
      if ( j > 0) yArraytr[j] += yArraytr[j-1];
    }

  for (int j=0;j<n;j++) yArraytr[j] /= yArraytr[n-1];
  ifFiletr.close();


}
//******************************************************
  /**
   * Destructor
   */
CMomDist::~CMomDist()
{
  delete [] xArrayz;
  delete [] yArrayz;
  delete [] xArraytr;
  delete [] yArraytr;
}
//*******************************************
  /**
   * Select Momentum components of the parent in a Monte Carlo fashion
   */
float CMomDist::getTransMom()
{
  float probtr = random.Rndm();
  int i = 0;
  for (;;)
    {
      if (yArraytr[i] > probtr) break;
      i++;
      if (i == n) break;
    }
  float xtr;
  if ( i==0) xtr = xArraytr[i];
  else xtr = xArraytr[i-1] + (xArraytr[i]-xArraytr[i-1])*random.Rndm();
  return xtr;
}

float CMomDist::getLongMom()
{
  float probz = random.Rndm();
  int i = 0;
  for (;;)
    {
      if (yArrayz[i] > probz) break;
      i++;
      if (i == n) break;
    }
  float xz;
  if ( i==0) xz = xArrayz[i];
  else xz = xArrayz[i-1] + (xArrayz[i]-xArrayz[i-1])*random.Rndm();
  return xz;
}
