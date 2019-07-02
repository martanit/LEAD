/*
 * integrator.h
 *
 *  Created on: March 11, 2019
 *  	Author: martanit
 */

#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

#include <cmath>
#include <memory>
#include <random>
#include "polymer.h"
#include "potential.h"
#include "extruder.h"
#include "vector_extruder.h"

class Integrator
{
  public: 
    Integrator(Polymer &poly, VectorExtruder &vector_extr, Parameters integrator_parm) : 
        m_vector_extr(vector_extr),
        m_integrator_timestep(integrator_parm.get_timestep()),
        m_integrator_gamma(integrator_parm.get_gamma()),
        m_integrator_temp(integrator_parm.get_temp()),
        gauss_term(0,1),
        transition_prob(0.,1.)
    {
      stoch_term =  std::sqrt(2.*m_integrator_gamma*
                                m_integrator_temp/
                                m_integrator_timestep);
      m_poly = std::make_unique<Polymer>(poly);
    };

    ~Integrator()
    {
    /*    delete m_poly.release();
        delete m_poly_new.release();
        delete m_poly_old.release();
        delete m_extr.release();
   */ };

    // Function to set and get polymer and vector extruder
    void set_new_polymer(Polymer& poly) { m_poly = std::make_unique<Polymer>(poly);}
    void set_new_extruder(VectorExtruder & new_vector_extr) { m_vector_extr = new_vector_extr; }
    const Polymer & get_poly() const {  return *m_poly; }
    const VectorExtruder & get_extr() const { return m_vector_extr; }
    
    // polymer integrators
    void langevin_overdamped();
    
    // extruder move
    void markov_chain();

  private:

    std::mt19937 mt {std::random_device{}()};    
    std::normal_distribution<double> gauss_term;
    std::uniform_real_distribution<double> transition_prob;
    
    std::unique_ptr<Polymer> m_poly;
    VectorExtruder m_vector_extr;
    
    double gauss_term_sphere;
    double stoch_term = 0;
    double m_integrator_timestep = 0.0001;
    double m_integrator_gamma=1.;
    double m_integrator_temp = 1.;
};
#endif /* INTEGRATOR_H_ */
