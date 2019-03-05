/*
 * Parameters.h
 *
 *  Created on: Feb 19, 2019
 *      Author: martanit
 */

#include <iostream> 	//std::cerr
#include <fstream>	//std::ifstream
#include <sstream>	//std::istringstream

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

class Parameters {

private:
	//dynamic parameters
	int m_nstep;		// number of step
	int m_print;		// verbosity
	float m_timestep;	// timestep

	//system parameters
	float m_temp;		// system temperature [K]
	std::string m_init;	// first configuration

	//polymer parameters
	float m_pmass;		// mass
	float m_pradius; 	// radius
	int m_psphere;		// number of interacting sphere
	float m_pdist;		// length of bond [angstrom]
	float m_bond;		// spring constant of harmonic oscillator (bond)
	float m_hradius;	//hard core radius of sphere

public:
	Parameters();
	Parameters( std::string );
	~Parameters();
	
	//function that read parameters from file
	bool read_parm( std::string );
	
	//function that store different type parameters
	template <typename myType>
	void set_parm( myType& m_parm, myType parm )  { m_parm = parm; }

	//function to access parameters
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
