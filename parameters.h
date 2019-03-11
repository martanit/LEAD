/*
 * parameters.h
 *
 *  Created on: Feb 19, 2019
 *      Author: martanit
 */

#include <iostream> 	// std::cerr
#include <fstream>	// std::ifstream
#include <sstream>	// std::istringstream

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

class Parameters {

private:
	// dynamic parameters
	int m_nstep = 100;		// number of step
	int m_print = 100;		// verbosity
	float m_timestep = 1.;		// timestep

	// system parameters
	float m_temp = 273.5;		// system temperature [K]
	std::string m_init = "init.dat";// first configuration

	// polymer parameters
	float m_pmass = 1.;		// mass
	float m_pradius = 0.1;	 	// radius
	int m_psphere = 100;		// number of interacting sphere
	float m_pdist = 10.;		// length of bond [angstrom]
	float m_bond = 15.;		// spring constant of harmonic oscillator (bond)
	float m_hradius = 5.;		// hard core radius of sphere

public:
	// default constructor
	Parameters() { };
	// constructor
	Parameters( std::string );
	~Parameters();
	
	// function that read parameters from file
	bool read_parm( std::string );
	
	// function that store different type parameters
	template <typename myType>
	void set_parm( myType& m_parm, myType parm )  { m_parm = parm; }

	// auxiliary function to access parameters
	int get_nstep() { return m_nstep; }
	int get_print() { return m_print; }
	float get_timestep() { return m_timestep; }

	float get_temp() { return m_temp; }
	
	float get_pmass() { return m_pmass; }
	int get_psphere() { return m_psphere; }
	float get_pdist() { return m_pdist; }
	float get_bond() { return m_bond; }
	float get_hradius() { return m_hradius; }
};

#endif /* PARAMETERS_H_ */
