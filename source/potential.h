/*
 * potential.h
 *
 *  Created on: March 11, 2019
 *  	Author: martanit
 */

#ifndef POTENTIAL_H_
#define POTENTIAL_H_

#include "parameters.h"
#include "polymer.h"
#include "vector_extruder.h"
#include <memory>

class Potential {

public:

    // construct potential from Polymer and set of parameters
    Potential(Polymer &poly, Parameters parm)
        : m_poly(poly),
	  sphere(poly.get_poly_nmonomers()),
          m_pot_epsilon(parm.get_epsilon()),
	  m_pot_rmin(parm.get_rmin()),
          m_pot_rcut(parm.get_rcut()),
	  m_pot_spring(poly.get_spring()),
          m_pot_box_length(parm.get_box_length()),
	  m_pot_kside(parm.get_kside()),
	  m_pot_eff_pot(parm.get_eff_pot()) {};

    Potential(Polymer &poly, VectorExtruder &vector_extr, Parameters parm)
        : m_poly(poly),
	  m_vector_extr(vector_extr),
	  sphere(poly.get_poly_nmonomers()),
	  extr_length(poly.get_poly_d()),
          m_pot_epsilon(parm.get_epsilon()),
	  m_pot_rmin(parm.get_rmin()),
          m_pot_rcut(parm.get_rcut()),
	  m_pot_spring(poly.get_spring()),
          m_pot_box_length(parm.get_box_length()),
	  m_pot_kside(parm.get_kside())	{};

    ~Potential() {};

    // access function
    void set_new_polymer(Polymer &new_poly) {
        m_poly = new_poly;
    }
    void set_new_extruder(VectorExtruder &new_vector_extr) {
        m_vector_extr = new_vector_extr;
    }

    const Polymer &get_poly() const {
        return m_poly;
    }
    const VectorExtruder &get_extr() const {
        return m_vector_extr;
    }

    // potentials
    void lennard_jones_f(int, bool);
    void well(bool);
    void soft_core_f(int, bool);
    void harmonic_spring_f(bool);
    void extruder_spring_f(bool);
    void box(bool);

protected:

    Polymer m_poly;
    VectorExtruder m_vector_extr;

private:


    // potential parameters
    double m_pot_epsilon = 1.;
    double m_pot_rmin = 1.;
    double m_pot_rmin_12 = std::pow(m_pot_rmin, 12);
    double m_pot_rmin_6 = std::pow(m_pot_rmin, 6);
    double m_pot_rcut = 5.;
    double m_pot_box_length = 10.;
    double m_pot_kside = 10.;
    double m_pot_spring = 250.;

    // set extruder bond as polymer bonds
    // with equilibrium distance as monomer
    // radius
    double k_extr = m_pot_spring;
    double extr_length = 0.89; 
    
    // potential interaction sphere
    std::vector<std::vector<int>> sphere;

    // effective potential from theory
    std::vector<double> m_pot_eff_pot;

    double dr = 0.;
    double x = 0., y = 0., z = 0.;
    double f_x = 0., f_y = 0., f_z = 0.;
    double spring_x = 0., spring_y = 0., spring_z = 0.;
    bool attractive = true;
};

#endif /* POTENTIAL_H_ */
