/*
 * parameters.h
 *
 *  Created on: Feb 19, 2019
 *      Author: martanit
 */

#include <fstream>  // std::ifstream
#include <iostream> // std::cerr
#include <sstream>  // std::istringstream
#include <vector>

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

class Parameters {
public:
  // default constructor
  Parameters(void){};
  // constructor
  Parameters(std::string, std::string);
  // constructor
  Parameters(std::string, std::string, std::string, std::string);
  ~Parameters();

  // function that read parameters, ctcf and
  // coupling probability (for extruder) from file
  bool read_parm(std::string);
  bool read_ctcf(std::string);
  bool print_ctcf(std::string);
  bool read_coupling_prob(std::string);
  bool print_param(std::string);

  // function that store different type parameters
  template <typename myType> void set_parm(myType &m_parm, myType parm) {
    m_parm = parm;
  }

  // auxiliary function to access parameters
  float get_nstep() { return m_nstep; }
  int get_print() { return m_print; }
  float get_timestep() { return m_timestep; }
  float get_gamma() { return m_gamma; }

  float get_temp() { return m_temp; }

  int get_nmonomers() { return m_nmonomers; }
  float get_diameter() { return m_diameter; }
  float get_spring() { return m_spring; }

  float get_epsilon() { return m_epsilon; }
  float get_rmin() { return m_rmin; }
  float get_rcut() { return m_rcut; }

  float get_permeability() { return m_perm_ctcf; }
  const std::vector<int> &get_ctcf() const { return m_ctcf; };
  const std::vector<double> &get_coupling_prob() const {
    return m_coupling_prob;
  };

  const float &get_rate_fwl() const { return m_rate_fwl; };
  const float &get_rate_fwr() const { return m_rate_fwr; };
  const float &get_rate_bwl() const { return m_rate_bwl; };
  const float &get_rate_bwr() const { return m_rate_bwr; };

  const float &get_kon() const { return m_k_on; }
  const float &get_koff() const { return m_k_off; }
  const float &get_k_diff() const { return m_k_diff; };
  const float &get_field_step() const { return m_field_step; };
  const float &get_field_length() const { return m_field_length; };
  const float &get_max_extr() const { return m_n_max_extr; }

protected:
  // dynamic parameters
  float m_nstep = 100;     // number of step
  int m_print = 100;     // verbosity
  float m_timestep = 1.; // timestep
  float m_gamma = 1.;    // parameter for langevin dynamics

  // system parameters
  float m_temp = 2.; // system temperature [K]

  // polymer parameters
  int m_nmonomers = 100;    // number of interacting monomers
  float m_diameter = 0.89; // diameter of monomers [50*nm]
  float m_spring = 500.;  // spring constant of harmonic oscillator (spring)
  std::string m_init = "init.dat"; // first configuration

  // potential parameters
  float m_epsilon = 1.; // parameter for LJ potential
  float m_rmin = 1.;    // parameter for LJ potential
  float m_rcut = 2.5;    // parameter for LJ potential action

  // extruder parameters
  float m_perm_ctcf = 0.1;
  std::vector<int> m_ctcf;
  std::vector<double> m_coupling_prob;

  float m_rate_fwl;
  float m_rate_fwr;
  float m_rate_bwl;
  float m_rate_bwr;

  // dynamics extruder parm
  float m_k_on = 0.5;
  float m_k_off = 0.001;

  // field of extruders parameters
  float m_k_diff = 1E3;
  float m_field_length = 10;
  float m_field_step = 10;
  float m_n_max_extr = 10;
};

#endif /* PARAMETERS_H_ */
