#include "sle.h"

/**
 *Constructor
 \param N0 is the number of equations
*/

sle::sle(int N0)
{
  N = N0;
  M = new float*[N];
  for (int i=0;i<N;i++) M[i] = new float[N];
  Y = new float[N];
}
//**************************************************************************
  /*
   * destructor
   */
sle::~sle()
{
  delete [] Y;
  for (int i=0;i<N;i++) delete [] M[i];
  delete [] M;
}
//**************************************************************************
  /**
   * Solves the equation, after one manually enters the values of M and Y
   */
void sle::solve()
{
  for (int j=0;j<N;j++)
    {
      float constant = M[j][j];
      for (int i=j;i<N;i++) M[j][i] /= constant;
      Y[j] /= constant;
      for (int i=0;i<N;i++)
	{
	  if (i != j)
	    {
              constant = M[i][j];
	      for (int k=j;k<N;k++) M[i][k] -= constant*M[j][k];
              Y[i] -= constant*Y[j];
	    }
	}
    }
}
