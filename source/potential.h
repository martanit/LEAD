/*
 * potential.h
 *
 *  Created on: March 11, 2019
 *  	Author: martanit
 */

#ifndef POTENTIAL_H_
#define POTENTIAL_H_

#include "polymer.h"
#include "extruder.h"

class Potential
{
  public:
    // default constructor
//    Potential()  { };
    
    // construct potential from Polymer and set of parameters
    Potential(Polymer& poly, std::vector<Extruder> extr, Parameters parm) : m_poly(poly),
                                                                m_extr(extr),
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
    void set_new_extruder(std::vector<Extruder> extr) { m_extr = extr;}
    const Polymer & get_poly() const {  return m_poly; }
    const std::vector<Extruder > get_extr() const {  return m_extr; }

    void lennard_jones_f();
    void harmonic_spring_f();
    void extruder_spring_f(); 

  protected:

    Polymer m_poly;
    std::vector<Extruder> m_extr;

    // potential parameters
    double m_pot_epsilon = 1.;
    double m_pot_sigma = 1.;
    double m_pot_rcut = 5.; 

    double k =  m_poly.get_bond(); 
    double k_extr = m_poly.get_bond()*5.;
    double extr_lenght = m_poly.get_poly_dist()/10.; 
    double m_box = 50.;
};

#endif /* POTENTIAL_H_ */
