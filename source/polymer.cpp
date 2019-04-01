#include "polymer.h"

Polymer::Polymer() : m_parm(Parameters())
{
}

Polymer::Polymer(std::string data) : m_parm(Parameters(data))
{
}

Polymer::Polymer
(std::vector<double>  polX, std::vector<double> polY, std::vector<double> polZ, std::string data):
  m_polX(polX), m_polY(polY), m_polZ(polZ), m_parm(Parameters(data))
{
}

Polymer::~Polymer()
{
}

void Polymer::first_sphere()
{	
  std::random_device rd;
  std::mt19937 mt (rd());
  std::uniform_real_distribution<double> dist(0., 1.);
	m_polX.push_back(dist(mt));
	m_polY.push_back(dist(mt));
	m_polZ.push_back(dist(mt));		
}

void Polymer::poly_configuration()
{
	double d = m_parm.get_pdist();
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> uniform01(0., 1.);
	this->first_sphere();
	for( int i = 1; i<=m_parm.get_psphere(); i++){
		bool is_overlap = true;
		while(is_overlap){
			
			double t = 2*M_PI*uniform01(mt);
			// uniform distribution on sphere
      double p = acos(1 - 2 * uniform01(mt));
				
			m_polX.push_back(m_polX[i-1]+d*std::sin(p)*std::cos(t)); 
			m_polY.push_back(m_polY[i-1]+d*std::sin(p)*std::sin(t));
			m_polZ.push_back(m_polZ[i-1]+d*std::cos(p));

			is_overlap = this->is_overlap(i);
		}
	}
}

bool Polymer::is_overlap(int lenght)
{
	bool is_overlap = false;
	for( int j = 0; j<lenght; j++){
		if( m_parm.get_hradius() > 
        this->dist(j, lenght)) 
      is_overlap = true;
		else is_overlap = false;		
	}
	return is_overlap;
}

bool Polymer::print_xyz( std::string out_xyz )
{	
  // open stream to write xyz file
  std::ofstream output;
  output.open( out_xyz, std::fstream::out );

	// return error if read file fail
	if( output.fail() ) {
		throw "ERROR: Impossible to write xyz file "+out_xyz;
		return 1;
	}
	
  output << m_parm.get_psphere() << std::endl;
	output << "Comment: prova"<< std::endl;
	for( int i = 0; i< m_parm.get_psphere(); i++){
		output << "Au "<< m_polX[i] << " " << m_polY[i] << " " << m_polZ[i] << std::endl;
	}
  output.close();
  return 0;
}

double Polymer::dist(int i, int j)
{
  return  sqrt(pow(m_polX[i]-m_polX[j], 2)+pow(m_polY[i]-m_polY[j], 2)+pow(m_polZ[i]-m_polZ[j], 2)); 
}
