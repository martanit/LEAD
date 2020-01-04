/*
 * integrator.h
 *
 *  Created on: March 11, 2019
 *  	Author: martanit
 */

#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

#include "extruder.h"
#include "polymer.h"
#include "potential.h"
#include "vector_extruder.h"
#include "cohesin_polymer.h"
#include <cmath>
#include <memory>
#include <random>

class Integrator {
public:
  Integrator(Polymer &poly, Parameters integrator_parm)
      : m_integrator_timestep(integrator_parm.get_timestep()),
        m_integrator_gamma(integrator_parm.get_gamma()),
        m_integrator_temp(integrator_parm.get_temp()), gauss_term(0, 1),
        transition_prob(0., 1.) {
    stoch_term = std::sqrt(2. * m_integrator_gamma * m_integrator_temp /
                           m_integrator_timestep);
    m_poly = std::make_unique<Polymer>(poly);
  };

  Integrator(Polymer &poly, VectorExtruder &vector_extr,
             Parameters integrator_parm)
      : m_vector_extr(vector_extr),
        m_integrator_timestep(integrator_parm.get_timestep()),
        m_integrator_gamma(integrator_parm.get_gamma()),
        m_integrator_temp(integrator_parm.get_temp()), gauss_term(0, 1),
        transition_prob(0., 1.) {
    stoch_term = std::sqrt(2. * m_integrator_gamma * m_integrator_temp /
                           m_integrator_timestep);
    m_poly = std::make_unique<Polymer>(poly);
  };

  Integrator(Polymer &poly, VectorExtruder &vector_extr,
	     CohesinPolymer & cohesin_field, Parameters integrator_parm)
      : m_vector_extr(vector_extr),
	m_cohesin_field(cohesin_field),
	new_cohesin_field(cohesin_field),
        m_integrator_timestep(integrator_parm.get_timestep()),
        m_integrator_gamma(integrator_parm.get_gamma()),
        m_integrator_temp(integrator_parm.get_temp()), gauss_term(0, 1),
        transition_prob(0., 1.), uniform05(0, 5) {
    stoch_term = std::sqrt(2. * m_integrator_gamma * m_integrator_temp /
                           m_integrator_timestep);
    m_poly = std::make_unique<Polymer>(poly);
  };

  ~Integrator(){};

  // Function to set and get polymer, vector extruder and field of extruders
  void set_new_polymer(Polymer &poly) {
    m_poly = std::make_unique<Polymer>(poly);
  }
  void set_new_extruder(VectorExtruder &new_vector_extr) {
    m_vector_extr = new_vector_extr;
  }
  void set_new_field(CohesinPolymer &new_cohesin_field){
    m_cohesin_field = new_cohesin_field;
  }

  const Polymer &get_poly() const { return *m_poly; }
  const VectorExtruder &get_extr() const { return m_vector_extr; }
  const CohesinPolymer &get_field() const { return m_cohesin_field; }

  // polymer integrators
  void langevin_overdamped();

  // extruder move
  void markov_chain();

  // field moves
  void extruders_diffusion();

private:
  // Random stuff
  std::mt19937 mt{std::random_device{}()};
  std::normal_distribution<double> gauss_term;
  std::uniform_real_distribution<double> transition_prob;
  std::uniform_int_distribution<int> uniform05;

  std::unique_ptr<Polymer> m_poly;
  VectorExtruder m_vector_extr;
  CohesinPolymer m_cohesin_field, new_cohesin_field;

  unsigned int direction;
  
  double gauss_term_sphere = 0;
  double stoch_term = 0;
  double m_integrator_timestep = 0.0001;
  double m_integrator_gamma = 1.;
  double m_integrator_temp = 1.;
};
#endif /* INTEGRATOR_H_ */
