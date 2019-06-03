/*
 * potential.h
 *
 *  Created on: March 11, 2019
 *  	Author: martanit
 */

#ifndef POTENTIAL_H_
#define POTENTIAL_H_

#include "polymer.h"
#include "vector_extruder.h"
#include<memory>

class Potential
{
  public:
    // construct potential from Polymer and set of parameters
    Potential(Polymer& poly, VectorExtruder& vector_extr, Parameters parm) : m_poly(poly),
                                                                m_vector_extr(vector_extr),
                                                                m_pot_epsilon(parm.get_epsilon()),
                                                                m_pot_sigma(parm.get_sigma()),
                                                                m_pot_rcut(parm.get_rcut()),
                                                                m_box(parm.get_box())
                                               
    {
    };
	 
    ~Potential() 
    { 
    };
   
    void set_new_polymer(Polymer& new_poly) { m_poly = new_poly;}
    void set_new_extruder(VectorExtruder& new_vector_extr) { m_vector_extr = new_vector_extr; }

    const Polymer & get_poly() const {  return m_poly; }
    const VectorExtruder & get_extr() const {  return m_vector_extr; }

    void lennard_jones_f();
    void harmonic_spring_f();
    void extruder_spring_f(); 

  protected:

    Polymer m_poly;
    VectorExtruder m_vector_extr;
    double x, y, z;
    unsigned n = 20000;
    unsigned r = n*(n-1)/2;

    std::vector<double> forces_local_x;
    std::vector<double> forces_local_y;
    std::vector<double> forces_local_z;
    
    // potential parameters
    double m_pot_epsilon = 1.;
    double m_pot_sigma = 1.;
    double m_pot_sigma_12 = std::pow(m_pot_sigma,12);
    double m_pot_sigma_6 = std::pow(m_pot_sigma,6);
    double m_pot_rcut = 5.; 

    double k =  m_poly.get_bond(); 
    double k_extr = m_poly.get_bond()*5.;
    double extr_lenght = m_poly.get_poly_dist()/10.; 
    double m_box = 50.;
 
    double dr;
    double spring_x = 0.,
           spring_y = 0.,
           spring_z = 0.;

    double lj_x = 0.,
           lj_y = 0.,
           lj_z = 0.;
};

#endif /* POTENTIAL_H_ */
