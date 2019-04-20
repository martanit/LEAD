#include "potential.h"

std::vector<double> Potential::lennard_jones(int i)
{
  std::vector<double> potential_vect;

  for( int j=0; j<m_poly.m_parm.get_psphere(); j++ ){
    if( j==i ) potential_vect.push_back(0) ;
    else potential_vect.push_back(
                 A/std::pow(m_poly.dist(j, i), 12)+ 
                 B/std::pow(m_poly.dist(j, i), 6));
  }
  return potential_vect;
}

std::vector<double> Potential::harmonic_spring()
{
  std::vector<double> bond_vect;
  for( int j=0; j < m_poly.m_parm.get_psphere(); j++ ){
    bond_vect.push_back( k * m_poly.m_parm.get_pdist() );
  }
  return bond_vect;
}
