#ifndef _sle
#define _sle


/**
 *!\brief solves a system of linear equations
 */
class sle
{
public:
  float ** M; //!< points to a matrix containing the coefficents of the equations
  float * Y; //!< points to a vector containing the solution 
  int N; //!< number of equations
  sle(int);
  ~sle();
  void solve();
};

#endif
