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
    // default constructor
    Potential()  { };
    
    // construct potential from Polymer and set of parameters
    Potential(Polymer& poly, Potential_Parameters parm) : m_poly(poly),
                                               m_conf(parm),
                                               m_pot_epsilon(parm.get_epsilon()),
                                               m_pot_sigma(parm.get_sigma()),
                                               m_pot_box(parm.get_box()),
                                               m_pot_rcut(parm.get_rcut())
                                               
    {
    };
	 
    ~Potential() 
    { 
    };
   
    void set_new_polymer(Polymer& poly) { m_poly = poly;}
    const Polymer & get_poly() const {  return m_poly; }

    void lennard_jones_f();
    void harmonic_spring_f();
      
  private:
    Polymer m_poly;
    Utils m_conf;
    
    // potential parameters
    double m_pot_epsilon = 1.;
    double m_pot_sigma = 1.;
    double m_pot_box = 10.;
    double m_pot_rcut = 5.; 

    double k =  m_poly.get_bond(); 
};

#endif /* POTENTIAL_H_ */
