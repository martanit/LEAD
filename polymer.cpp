/*
 * Polymer.cpp
 *
 *  Created on: Feb 19, 2019
 *      Author: martanit
 */

#include "polymer.h"



Polymer::Polymer(std::array<double, int psphere> polX, std::array<double, int psphere> polY, std::array<double, int psphere> polZ) : m_polX(polX), m_polY(polY), m_polZ(polZ) {
	m_psphere = psphere;
	m_pdist = 	
}

Polymer::~Polymer() {
	// TODO Auto-generated destructor stub
}

void Polymer::first_sphere() {
	
 	std::random_device rd;
    	std::mt19937 mt(rd());
    	std::uniform_real_distribution<double> dist(0., 1.);

	m_polX[0] = dist(mt);
	m_polY[0] = dist(mt);
	m_polz[0] = dist(mt);
		
}

void Polymer::poly_configuration() {
	
	d = m_parm.get_pdist();

	first_sphere();
	for( int i = 1; i<m_parm.get_psphere(); i++){		
		bool is_overlap = true;
		while(is_overlap){
 			std::random_device rd;
    			std::mt19937 mt(rd());
    			std::uniform_real_distribution<double> theta(0., 2.*M_PI);
			std::uniform_real_distribution<double> phi(0., M_PI);

			m_polX[i] = m_polX[i-1]+d*sin(phi)*cos(theta); 
			m_polY[i] = m_polY[i-1]+d*sin(phi)*sin(theta);
			m_polZ[i] = m_polZ[i-1]+d*cos(phi);

			is_overlap = is_overlap(&i);
		}
	}
}

bool Polymer::is_overlap(int & i){
	bool is_overlap = false;
	r_h = m_parm.get_hardcore()
	for( int j = 1; j<i; j++){
		if( m_parm.get_hardcore() > sqrt(pow(m_polX[j]-m_polX[i], 2)+pow(m_polY[j]-m_polY[i], 2)+pow(m_polZ[j]-m_polZ[i], 2))) is_overlap = true;
		else is_overlap = false;		
	}
	return is_overlap;
}	

