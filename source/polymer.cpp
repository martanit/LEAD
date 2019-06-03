#include "polymer.h"

bool print_xyz(Polymer&, std::string);
bool read_xyz(Polymer&, std::string);

Polymer::Polymer()
{
    this->set_size(); 
    this->reset_force();
    this->poly_configuration();
    this->poly_velocity();
}

Polymer::Polymer(Parameters parm) : m_poly_mass(parm.get_pmass()),
                                    m_poly_sphere(parm.get_psphere()),
                                    m_poly_dist(parm.get_pdist()),
                                    m_poly_bond(parm.get_bond()),
                                    m_poly_hradius(parm.get_hradius()),
                                    m_box(parm.get_box()),
                                    uniform01(0.,1.),
                                    uniform0505(-0.5, 0.5)
{
  this->set_size(); 
  this->reset_force();
  this->poly_configuration();
  this->poly_velocity();
}

Polymer::Polymer (Parameters parm, std::string poly_xyz) : m_poly_mass(parm.get_pmass()),
                                                           m_poly_sphere(parm.get_psphere()),
                                                           m_poly_dist(parm.get_pdist()),
                                                           m_poly_bond(parm.get_bond()),
                                                           m_poly_hradius(parm.get_hradius()),
                                                           m_box(parm.get_box()),
                                                           uniform01(0.,1.),
                                                           uniform0505(-0.5, 0.5)
{
   this->set_size(); 
   this->reset_force();
   read_xyz((*this), poly_xyz);
  // this->control_poly;
}

Polymer::~Polymer()
{
}

void Polymer::set_size()
{
    m_poly_x.resize(m_poly_sphere);  
    m_poly_y.resize(m_poly_sphere); 
    m_poly_z.resize(m_poly_sphere); 
    
    m_poly_vx.resize(m_poly_sphere);
    m_poly_vy.resize(m_poly_sphere); 
    m_poly_vz.resize(m_poly_sphere); 
    
    m_poly_fx.resize(m_poly_sphere); 
    m_poly_fy.resize(m_poly_sphere); 
    m_poly_fz.resize(m_poly_sphere);
}

void Polymer::first_sphere()
{	
	m_poly_x[0] = pbc( uniform01(mt) );   // x
	m_poly_y[0] = pbc( uniform01(mt) );   // y
	m_poly_z[0] = pbc( uniform01(mt) );	// z	
}

void Polymer::poly_configuration()
{
    double d = m_poly_dist;
	this->first_sphere();
	for(unsigned int i = 1; i<m_poly_sphere; i++){
		bool is_overlap = true;
		while(is_overlap){
			
			double t = 2*M_PI*uniform01(mt);
			// uniform distribution on sphere
            double p = acos(1 - 2 * uniform01(mt));
      
            m_poly_x[i] = pbc( m_poly_x[i-1]+d*std::sin(p)*std::cos(t) );
            m_poly_y[i] = pbc( m_poly_y[i-1]+d*std::sin(p)*std::sin(t) );
            m_poly_z[i] = pbc( m_poly_z[i-1]+d*std::cos(p) );
			is_overlap = this->is_overlap(i);
		}
	}
}

bool Polymer::is_overlap(int lenght)
{
    bool is_overlap = false;
	for(unsigned int j = 0; j<lenght; j++){
		if( this->dist(j, lenght) < (2*m_poly_hradius)) {
            is_overlap = true;
            break;
        }
		else 
            is_overlap = false;		
        }
	return is_overlap;
}

void Polymer::poly_velocity()
{
    std::vector<double> sum(3, 0.0);

    for( unsigned int i = 0; i<m_poly_sphere; i++){
        m_poly_vx[i] =  uniform0505(mt);
        m_poly_vy[i] =  uniform0505(mt);
        m_poly_vz[i] =  uniform0505(mt);
      
        sum[0] += m_poly_vx[i];    
        sum[1] += m_poly_vy[i];    
        sum[2] += m_poly_vz[i];    
    }
  
    for( unsigned short int i=0; i<3; i++) 
        sum[i] /= m_poly_sphere;
   
    for( unsigned int i = 0; i<m_poly_sphere; i++){
        m_poly_vx[i] -= sum[0];
        m_poly_vy[i] -= sum[1];
        m_poly_vz[i] -= sum[2];
    }
}

double Polymer::dist(const int& i, const int& j)
{
    
    d = std::sqrt(std::pow(m_poly_x[i]-m_poly_x[j], 2)+
                            std::pow(m_poly_y[i]-m_poly_y[j], 2)+
                            std::pow(m_poly_z[i]-m_poly_z[j], 2));
    return (d < 10E-9) ? 10E-9 : d;
}

void Polymer::reset_force()
{
    std::fill(m_poly_fx.begin(), m_poly_fx.end(), 0.);
    std::fill(m_poly_fy.begin(), m_poly_fy.end(), 0.);
    std::fill(m_poly_fz.begin(), m_poly_fz.end(), 0.);
}

void Polymer::add_force(const int& i, const double& fx, const double& fy, const double& fz)
{
    m_poly_fx[i] += fx;
    m_poly_fy[i] += fy;
    m_poly_fz[i] += fz;
}

void Polymer::reset_energy()
{
   m_poly_e=0.;
}

void Polymer::add_energy(const double& e)
{
    m_poly_e += e;
}
// Input/Output function for read and write polymers
bool print_xyz(Polymer& poly, std::string out_xyz )
{	
  // open stream to write xyz file
  std::ofstream output;
  output.open( out_xyz, std::fstream::app );

	// return error if read file fail
	if( output.fail() ) {
		throw "ERROR: Impossible to write xyz file "+out_xyz;
		return 1;
	}
	
  output << poly.m_poly_sphere << std::endl << std::endl;
	for( int i = 0; i< poly.m_poly_sphere; i++){
        output << "\tAu\t\t\t" << poly.m_poly_x[i] << "\t\t\t" << poly.m_poly_y[i] << "\t\t\t" << poly.m_poly_z[i]  << std::endl;
	}

  output.close();
  return 0;
}


bool read_xyz(Polymer& poly, std::string in_xyz)
{
  std::string name;
  std::ifstream read_conf;
  read_conf.open( in_xyz, std::fstream::in );
  if( read_conf.fail() ){
		throw "ERROR: Impossible to open polymer coordinates file "+in_xyz;
        return 1;
    }
  
  std::string dummyLine;
  getline(read_conf, dummyLine);
  getline(read_conf, dummyLine);
  
  std::string line;
  int i=0;
   while(std::getline(read_conf, line)){
    std::stringstream iss(line);
    iss >> name >> poly.m_poly_x[i]  >> poly.m_poly_y[i]  >> poly.m_poly_z[i] >> poly.m_poly_vx[i] >> poly.m_poly_vy[i] >> poly.m_poly_vz[i];
    i++;
   }
  read_conf.close();
  return 0;
}
