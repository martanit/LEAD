/*
 * parameters.cpp
 *
 *  Created on: Feb 19, 2019
 *      Author: martanit
 */

#include "parameters.h"

Parameters::Parameters( std::string file_name, std::string ctcf, std::string prob, std::string rate )
{
	this -> read_parm( file_name );
    this -> read_ctcf( ctcf );
    this -> read_coupling_prob( prob );
    this -> read_rate( rate );
}

Parameters::~Parameters()
{
}

bool Parameters::read_parm( std::string file_name )
{
	// read file
	std::ifstream parm_file;
	parm_file.open( file_name, std::fstream::in );

	// return error if read file fail
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
                else if( name == "box" ) set_parm(m_box, std::stof(par));
				else if( name == "init" ) set_parm(m_init, par);
				else if( name == "psphere") set_parm( m_psphere, std::stoi(par));
				else if( name == "pmass" ) set_parm(m_pmass, std::stof(par));
				else if( name == "pdist" ) set_parm(m_pdist, std::stof(par));		
				else if( name == "bond" ) set_parm(m_bond, std::stof(par));
				else if( name == "hradius" ) set_parm(m_hradius, std::stof(par));
				else if( name == "epsilon" ) set_parm(m_epsilon, std::stof(par));
				else if( name == "sigma" ) set_parm(m_sigma, std::stof(par));
                else if( name == "rcut" ) set_parm(m_rcut, std::stof(par));
                else if( name == "gamma" ) set_parm(m_gamma, std::stof(par));
                else if( name == "k_on" ) set_parm(m_k_on, std::stof(par));
                else if( name == "k_off" ) set_parm(m_k_off, std::stof(par));
                else if( name == "n_max_extr" ) set_parm(m_n_max_extr, std::stoi(par));
			}
		}
	}
	parm_file.close();
	return 0;
}

bool Parameters::read_ctcf( std::string file_name)
{
    int x;

	// read file
	std::ifstream ctcf_file;
	ctcf_file.open( file_name, std::fstream::in );

	// return error if read file fail
	if( ctcf_file.fail() ) {
		throw "ERROR: Impossible to open parameters file "+file_name;
		return 1;
	}

	std::string line;
	while ( std::getline(ctcf_file, line) ) {
		std::istringstream iss(line);
	    iss >> x;		
        m_ctcf.push_back(x);
    }
	ctcf_file.close();
	return 0;
}

bool Parameters::read_coupling_prob( std::string file_name)
{
    double p;

	// read file
	std::ifstream coupling_prob_file;
	coupling_prob_file.open( file_name, std::fstream::in );

	// return error if read file fail
	if( coupling_prob_file.fail() ) {
		throw "ERROR: Impossible to open parameters file "+file_name;
		return 1;
	}

	std::string line;
	while ( std::getline(coupling_prob_file, line) ) {
		std::istringstream iss(line);
	    iss >> p;		
        m_coupling_prob.push_back(p);
    }
	coupling_prob_file.close();
	return 0;
}

bool Parameters::read_rate( std::string file_name)
{
    double l;
    double r;

	// read file
	std::ifstream rate_file;
	rate_file.open( file_name, std::fstream::in );

	// return error if read file fail
	if( rate_file.fail() ) {
		throw "ERROR: Impossible to open parameters file "+file_name;
		return 1;
	}

	std::string line;
	while ( std::getline(rate_file, line) ) {
		std::istringstream iss(line);
	    iss >> l >> r;		
        m_rate_l.push_back(l);
        m_rate_r.push_back(r);
    }
	rate_file.close();
	return 0;
}
