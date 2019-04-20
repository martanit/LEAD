#include "polymer.h"

Polymer::Polymer() :
                m_poly_r_v(m_poly_mass, 6)
{
  this->poly_configuration();
}

Polymer::Polymer(Parameters parm) : 
                      m_poly_r_v(parm.get_psphere(),6),
                      m_poly_mass(parm.get_pmass()),
                      m_poly_sphere(parm.get_psphere()),
                      m_poly_dist(parm.get_pdist()),
                      m_poly_bond(parm.get_bond()),
                      m_poly_hradius(parm.get_hradius())
{
  this->poly_configuration();
}

Polymer::Polymer (Parameters parm, std::string poly_xyz) :
                      m_poly_r_v(parm.get_psphere(),6),
                      m_poly_mass(parm.get_pmass()),
                      m_poly_sphere(parm.get_psphere()),
                      m_poly_dist(parm.get_pdist()),
                      m_poly_bond(parm.get_bond()),
                      m_poly_hradius(parm.get_hradius())
                      
{
  this->read_xyz(poly_xyz);
  //this->control_poly;
}

Polymer::~Polymer()
{
}

void Polymer::first_sphere()
{	
  std::random_device rd;
  std::mt19937 mt (rd());
  std::uniform_real_distribution<double> dist(0., 1.);
	m_poly_r_v(0,0) = pbc( dist(mt) ); // x
	m_poly_r_v(0,1) = pbc( dist(mt) ); // y
	m_poly_r_v(0,2) = pbc( dist(mt) );	// z	
}

void Polymer::poly_configuration()
{
	double d = m_poly_dist;
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> uniform01(0., 1.);
	this->first_sphere();
	for( int i = 1; i<m_poly_sphere; i++){
		bool is_overlap = true;
		while(is_overlap){
			
			double t = 2*M_PI*uniform01(mt);
			// uniform distribution on sphere
      double p = acos(1 - 2 * uniform01(mt));
				
      m_poly_r_v(i,0) = pbc( m_poly_r_v(i-1,0)+d*std::sin(p)*std::cos(t) );
      m_poly_r_v(i,1) = pbc( m_poly_r_v(i-1,1)+d*std::sin(p)*std::sin(t) );
      m_poly_r_v(i,2) = pbc( m_poly_r_v(i-1,2)+d*std::cos(p) );
			
			is_overlap = this->is_overlap(i);
		}
	}
}

bool Polymer::is_overlap(int lenght)
{
	bool is_overlap = false;
	for( int j = 0; j<lenght; j++){
		if( m_poly_hradius > 
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
	
  output << m_poly_sphere << std::endl;
	output << "Comment: prova"<< std::endl;
	for( int i = 0; i< m_poly_sphere; i++){
		output << "Au "<< m_poly_r_v(i,0) << " " << m_poly_r_v(i,1) << " " << m_poly_r_v(i,2) << std::endl;
	}
  output.close();
  return 0;
}

bool Polymer::read_xyz( std::string in_xyz)
{

  std::ifstream read_conf;
  read_conf.open( in_xyz, std::fstream::in );
  if( read_conf.fail() ){
		throw "ERROR: Impossible to open polymer coordinates file "+in_xyz;
		return 1;
	}
  std::string line;
  read_conf.seekg(2);
  int i=0;
  while(std::getline(read_conf, line)){
    std::stringstream iss(line);
    iss  >> m_poly_r_v(i,0) >> m_poly_r_v(i,1) >> m_poly_r_v(i, 2);   
    i++;
   }
  read_conf.close();
  return 0;
}

double Polymer::dist(int i, int j)
{
  return  std::sqrt(std::pow(m_poly_r_v(i,0)-m_poly_r_v(j,0), 2)+
                    std::pow(m_poly_r_v(i,1)-m_poly_r_v(j,1), 2)+
                    std::pow(m_poly_r_v(i,2)-m_poly_r_v(j,2), 2)); 
}

//Algorithm for periodic boundary conditions with side L=box
double Polymer::pbc(double r)
{
    return r - 10. * rint(r/10.);
}
