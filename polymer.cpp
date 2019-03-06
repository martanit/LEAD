#include "polymer.h"


Polymer::Polymer(){};
Polymer::Polymer(std::vector<double>  polX, std::vector<double> polY, std::vector<double> polZ) : m_polX(polX), m_polY(polY), m_polZ(polZ) {
	m_psphere = m_parm.get_psphere();
}

Polymer::~Polymer() {
}

void Polymer::first_sphere() {
	
 	std::random_device rd;
    	std::mt19937 mt(rd());
    	std::uniform_real_distribution<double> dist(0., 1.);

	m_polX.push_back(dist(mt));
	m_polY.push_back(dist(mt));
	m_polZ.push_back(dist(mt));
		
}

void Polymer::poly_configuration() {
	double d = m_parm.get_pdist();
	this->first_sphere();
	std::cout << m_parm.get_psphere() << std::endl;
	std::cout<< "Comment: prova"<< std::endl;
	for( int i = 1; i<=m_parm.get_psphere(); i++){		
		bool is_overlap = true;
		while(is_overlap){
 			std::random_device rd;
    			std::mt19937 mt(rd());
    			std::uniform_real_distribution<double> theta(0., 2.*M_PI);
			std::uniform_real_distribution<double> phi(0., M_PI);
			
			double t = theta(mt);
			double p = phi(mt);
				
			m_polX.push_back(m_polX[i-1]+d*std::sin(p)*std::cos(t)); 
			m_polY.push_back(m_polY[i-1]+d*std::sin(p)*std::sin(t));
			m_polZ.push_back(m_polZ[i-1]+d*std::cos(p));

			is_overlap = this->is_overlap(i);
		}
		std::cout << "Au "<< m_polX[i] << " " << m_polY[i] << " " << m_polZ[i] << std::endl;
	}
}

bool Polymer::is_overlap(int lenght){
	bool is_overlap = false;
	for( int j = 0; j<lenght; j++){
		if( m_parm.get_hradius() > sqrt(pow(m_polX[j]-m_polX[lenght], 2)+pow(m_polY[j]-m_polY[lenght], 2)+pow(m_polZ[j]-m_polZ[lenght], 2))) is_overlap = true;
		else is_overlap = false;		
	}
	return is_overlap;
}	

