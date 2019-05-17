/*
 * potential.h
 *
 *  Created on: March 11, 2019
 *  	Author: martanit
 */

#ifndef POTENTIAL_H_
#define POTENTIAL_H_

#include "polymer.h"
#include "utils.h"

class Potential
{
  public:
    // default constructor
    Potential()  { };
    
    // construct potential from Polymer and set of parameters
    Potential(Polymer& poly, Parameters parm) : m_poly(poly),
                                                m_pot_epsilon(parm.get_epsilon()),
                                                m_pot_sigma(parm.get_sigma()),
                                                m_pot_rcut(parm.get_rcut()),
                                                m_box(parm.get_box())
                                               
    {
    };
	 
    ~Potential() 
    { 
    };
   
    void set_new_polymer(Polymer& poly) { m_poly = poly;}
    const Polymer & get_poly() const {  return m_poly; }

    void lennard_jones_f();
    void harmonic_spring_f();

    friend const double pbc(double r);

  private:
    Polymer m_poly;
    
    // potential parameters
    double m_pot_epsilon = 1.;
    double m_pot_sigma = 1.;
    double m_pot_rcut = 5.; 

    double k =  m_poly.get_bond(); 
    double m_box = 50.;
};

#endif /* POTENTIAL_H_ */
