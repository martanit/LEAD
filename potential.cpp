#include "potential.h"


std::vector<double> Potential::lennard_jones(int i)
{
  std::vector<double> potential_vect;

  for( int j=0; j<m_poly.m_parm.get_psphere(); j++){
     potential_vect.push_back(
         4 * m_epsilon *
            (pow((m_poly.m_parm.get_pdist()/
                  m_poly.dist(j, i)), 12) 
           - pow((m_poly.m_parm.get_pdist()/
                  m_poly.dist(j, i)), 6)));
  }
  return potential_vect;
}

double HarmonicSpring(double d)
{
  return 1/2 * k * d;
}
