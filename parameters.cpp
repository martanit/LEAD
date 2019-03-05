/*
 * Parameters.cpp
 *
 *  Created on: Feb 19, 2019
 *      Author: martanit
 */

#include "parameters.h"

//default constructor
Parameters::Parameters() {
	m_nstep = 100;
	m_print = 100;
	m_timestep = 1.;
	m_temp = 273.5;
	m_init = "init.dat";

	m_psphere = 100;
	m_pradius = 1.;
	m_pmass = 1.;
	m_pdist = 1.;
	m_bond = 1.;
	m_hradius = 1.;
}

//constructor
Parameters::Parameters( std::string file_name ){
	this -> read_parm( file_name );
}

Parameters::~Parameters() { }

bool Parameters::read_parm( std::string file_name ) {

	//read file
	std::ifstream parm_file;
	parm_file.open( file_name, std::fstream::in );

	//return error if read file fail
	if( parm_file.fail() ) {
		throw "ERROR: Impossible to open parameters file "+file_name;
		return 1;
	}

	std::string line;
	while ( std::getline(parm_file, line) ) {
		std::istringstream is_line(line);
		std::string name;
		if ( std::getline(is_line, name, '=') ) {
			std::string par;
			if ( std::getline( is_line, par ) ) {
				if( name == "nstep" ) set_parm(m_nstep, std::stoi(par));
				else if( name == "print" ) set_parm(m_print, std::stoi(par));
				else if( name == "timestep" ) set_parm(m_timestep, std::stof(par));	
				else if( name == "temp" ) set_parm(m_temp, std::stof(par));
				else if( name == "init" ) set_parm(m_init, par);
				else if( name == "m_psphere") set_parm( m_psphere, std::stoi(par));
				else if( name == "m_pradius") set_parm( m_pradius, std::stof(par));
				else if( name == "pmass" ) set_parm(m_pmass, std::stof(par));
				else if( name == "pdist" ) set_parm(m_pdist, std::stof(par));		
				else if( name == "bond" ) set_parm(m_bond, std::stof(par));
				else if( name == "hradius" ) set_parm(m_hradius, stof(par));
			}
		}
	}
	parm_file.close();
	return 0;
}
