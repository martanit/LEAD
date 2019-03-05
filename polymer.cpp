#include "polymer.h"

Polymer::Polymer(std::vector<double>  polX, std::vector<double> polY, std::vector<double> polZ) : m_polX(polX), m_polY(polY), m_polZ(polZ) {
	m_psphere = m_parm.get_psphere();
}

Polymer::~Polymer() {
}

void Polymer::first_sphere() {
	
 	std::random_device rd;
    	std::mt19937 mt(rd());
    	std::uniform_real_distribution<double> dist(0., 1.);

	m_polX[0] = dist(mt);
	m_polY[0] = dist(mt);
	m_polZ[0] = dist(mt);
		
}

void Polymer::poly_configuration() {
	
	double d = m_parm.get_pdist();
	first_sphere();
	for( int i = 1; i<m_parm.get_psphere(); i++){		
		bool is_overlap = true;
		while(is_overlap){
 			std::random_device rd;
    			std::mt19937 mt(rd());
    			std::uniform_real_distribution<double> theta(0., 2.*M_PI);
			std::uniform_real_distribution<double> phi(0., M_PI);
			
			m_polX[i] = m_polX[i-1]+d*std::sin(phi(mt))*std::cos(theta(mt)); 
			m_polY[i] = m_polY[i-1]+d*std::sin(phi(mt))*std::sin(theta(mt));
			m_polZ[i] = m_polZ[i-1]+d*std::cos(phi(mt));
			is_overlap = this->is_overlap(i);
		}
	}
}

bool Polymer::is_overlap(int lenght){
	bool is_overlap = false;
	for( int j = 1; j<lenght; j++){
		if( m_parm.get_hradius() > sqrt(pow(m_polX[j]-m_polX[lenght], 2)+pow(m_polY[j]-m_polY[lenght], 2)+pow(m_polZ[lenght]-m_polZ[lenght], 2))) is_overlap = true;
		else is_overlap = false;		
	}
	return is_overlap;
}	

