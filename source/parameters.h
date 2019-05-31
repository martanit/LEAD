/*
 * parameters.h
 *
 *  Created on: Feb 19, 2019
 *      Author: martanit
 */

#include <iostream> 	// std::cerr
#include <fstream>	// std::ifstream
#include <sstream>	// std::istringstream
#include <vector>

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

class Parameters
{
public:
	// default constructor
	Parameters(void) { };
	// constructor
	Parameters( std::string, std::string, std::string, std::string );
	Parameters( std::string, std::string, std::string );
	~Parameters();
	
	// function that read parameters, ctcf and 
    // coupling probability (for extruder) from file
	bool read_parm( std::string );
    bool read_ctcf( std::string );
    bool read_coupling_prob( std::string );
    bool read_rate( std::string );

	// function that store different type parameters
	template <typename myType>
	void set_parm( myType& m_parm, myType parm )  { m_parm = parm; }

	// auxiliary function to access parameters
	int get_nstep() { return m_nstep; }
	int get_print() { return m_print; }
    float get_timestep() { return m_timestep; }
    float get_gamma() { return m_gamma; }
    
    float get_temp() { return m_temp; }
	float get_box() { return m_box; }
    
    float get_pmass() { return m_pmass; }
	int get_psphere() { return m_psphere; }
	float get_pdist() { return m_pdist; }
	float get_bond() { return m_bond; }
	float get_hradius() { return m_hradius; }
    
    float get_epsilon() { return m_epsilon; }
    float get_sigma() { return m_sigma; }
    float get_rcut() { return m_rcut; }
    
    float get_permeability() { return m_perm_ctcf; } 
    const std::vector<int>& get_ctcf() const { return m_ctcf; };
    const std::vector<double>& get_coupling_prob() const { return m_coupling_prob; };
   
    const float& get_rate_l() const { return m_rate_l; };
    const float& get_rate_r() const { return m_rate_r; };

    const float& get_kon() const {return m_k_on; }
    const float& get_koff() const {return m_k_off; }
    const int& get_max_extr() const {return m_n_max_extr;}

    const std::vector<double>& get_rate_vl() { return m_rate_vl; };
    const std::vector<double>& get_rate_vr() { return m_rate_vr; };

protected:
	// dynamic parameters
	int m_nstep = 100;		// number of step
	int m_print = 100;		// verbosity
	float m_timestep = 1.;		// timestep
    float m_gamma = 1.; // parameter for langevin dynamics

	// system parameters
	float m_temp = 2.;		// system temperature [K]
    float m_box = 10.;

	// polymer parameters
	float m_pmass = 1.;		// mass
	int m_psphere = 10;		// number of interacting sphere
	float m_pdist = 10.;		// length of bond [angstrom]
	float m_bond = 15.;		// spring constant of harmonic oscillator (bond)
	float m_hradius = 5.;		// hard core radius of sphere
	std::string m_init = "init.dat";// first configuration
 
    // potential parameters
    float m_epsilon = 1.; // parameter for LJ potential
    float m_sigma = 1.; // parameter for LJ potential
    float m_rcut = 5.; // parameter for LJ potential action
    
    // extruder parameters
    float m_perm_ctcf = 0.9;
    std::vector<int> m_ctcf;
    std::vector<double> m_coupling_prob;
    float m_rate_l;
    float m_rate_r;
    
    //use for different rate for each sphere
    std::vector<double> m_rate_vl;
    std::vector<double> m_rate_vr;

    // dynamics extruder parm
    float m_k_on = 0.5;
    float m_k_off = 0.998;
    int m_n_max_extr = 10;
};

#endif /* PARAMETERS_H_ */
