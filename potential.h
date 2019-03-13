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
	  Potential(Polymer & poly) : m_poly(poly) { };
	  ~Potential() { };
    
    // function that return lennard jones potential 
    // between i-atom and others
    std::vector<double> lennard_jones(int );
	
    double HarmonicSpring(double );

  protected:
	  Polymer m_poly;
    double m_epsilon;
};

#endif /* POTENTIAL_H_ */
