#include "fragment.h"

/**
 * Constructor
  \param Z0 is the atomic number of the fragment
  \param mass0 is the mass of teh fragment in amu
 */
fragment::fragment(float Z0,float mass0)
{
  Z = Z0;
  mass = mass0;
}
/**
 *Initialize a new fragment
  \param Z0 is the atomic number of the fragment
  \param mass0 is the mass of teh fragment in amu
 */
void fragment::init(float Z0,float mass0)
{
  Z = Z0;
  mass = mass0;
}
//*********************************************************
  /**
   * Destructor
   */
fragment::~fragment()
{

}

