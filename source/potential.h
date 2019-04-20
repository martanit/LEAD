/*
 * potential.h
 *
 *  Created on: March 11, 2019
 *  	Author: martanit
 */

#ifndef POTENTIAL_H_
#define POTENTIAL_H_

#include "polymer.h"

class Potential
{
  public:
	  Potential(Parameters parm, Polymer & poly) : m_poly(poly) { };
	  ~Potential() { };
    
    // function that return lennard jones potential 
    std::vector<double> lennard_jones(int );

    // spring chain forces  
    std::vector<double> harmonic_spring( );

  private:
	  Polymer m_poly;

    // function constant
    const double A = 4*m_poly.m_parm.get_epsilon()*std::pow(m_poly.m_parm.get_sigma(),12);
    const double B = 4*m_poly.m_parm.get_epsilon()*std::pow(m_poly.m_parm.get_sigma(),6);
    const double k =  1./2. * m_poly.m_parm.get_bond(); 
};

#endif /* POTENTIAL_H_ */
