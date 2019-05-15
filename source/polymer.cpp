#include "polymer.h"

Polymer::Polymer() :
                m_poly_r_v(m_poly_sphere=10, 6),
                m_conf()
{
  m_poly_force.zeros();
  this->poly_configuration();
  this->poly_velocity();
}

Polymer::Polymer(Parameters parm) : 
                      m_poly_r_v(parm.get_psphere(),6),
                      m_poly_force(parm.get_psphere(),3),
                      m_conf(parm),
                      m_poly_mass(parm.get_pmass()),
                      m_poly_sphere(parm.get_psphere()),
                      m_poly_dist(parm.get_pdist()),
                      m_poly_bond(parm.get_bond()),
                      m_poly_hradius(parm.get_hradius())
{
  m_poly_force.zeros();
  this->poly_configuration();
  this->poly_velocity();
}

Polymer::Polymer (Parameters parm, std::string poly_xyz) :
                      m_poly_r_v(parm.get_psphere(),6),
                      m_poly_force(parm.get_psphere(),3),
                      m_conf(parm),
                      m_poly_mass(parm.get_pmass()),
                      m_poly_sphere(parm.get_psphere()),
                      m_poly_dist(parm.get_pdist()),
                      m_poly_bond(parm.get_bond()),
                      m_poly_hradius(parm.get_hradius())
                      
{
   m_poly_force.zeros();
   m_poly_r_v = read_xyz(poly_xyz);
  // this->control_poly;
}

Polymer::~Polymer()
{
}

void Polymer::first_sphere()
{	
  std::random_device rd;
  std::mt19937 mt (rd());
  std::uniform_real_distribution<double> dist(0., 1.);
	m_poly_r_v(0,0) = m_conf.pbc( dist(mt) ); // x
	m_poly_r_v(0,1) = m_conf.pbc( dist(mt) ); // y
	m_poly_r_v(0,2) = m_conf.pbc( dist(mt) );	// z	
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
      
      m_poly_r_v(i,0) = m_conf.pbc( m_poly_r_v(i-1,0)+d*std::sin(p)*std::cos(t) );
      m_poly_r_v(i,1) = m_conf.pbc( m_poly_r_v(i-1,1)+d*std::sin(p)*std::sin(t) );
      m_poly_r_v(i,2) = m_conf.pbc( m_poly_r_v(i-1,2)+d*std::cos(p) );
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

void Polymer::poly_velocity()
{
  std::vector<double> sum(3, 0.0);

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> uniform01(-0.5,0.5);
  for( unsigned int i = 0; i<m_poly_sphere; i++){
    for( auto k : {3,4,5}) m_poly_r_v(i,k) = uniform01(mt);
    for( auto k : {3,4,5}) sum[k-3] +=  m_poly_r_v(i,k);
  }
  
  for( unsigned short int i=0; i<3; i++) 
    sum[i] /= m_poly_sphere;
   
  for( int i = 0; i<m_poly_sphere; i++){
    for( auto k : {3,4,5}) m_poly_r_v(i,k) = m_poly_r_v(i,k)-sum[k-3];
  }
}

double Polymer::dist(int i, int j)
{
    
    double dist = std::sqrt(std::pow(m_poly_r_v(i,0)-m_poly_r_v(j,0), 2)+
                            std::pow(m_poly_r_v(i,1)-m_poly_r_v(j,1), 2)+
                            std::pow(m_poly_r_v(i,2)-m_poly_r_v(j,2), 2));
    if(dist < 0.0000000001) return 0.0000000001;
    else return dist;
}

void Polymer::set_force(int i, double fx, double fy, double fz)
{
  m_poly_force(i, 0) = fx;
  m_poly_force(i, 1) = fy;
  m_poly_force(i, 2) = fz;
}

void Polymer::reset_force()
{
  m_poly_force.zeros();
}

void Polymer::add_force(int i, double f)
{
  m_poly_force(i, 0) += f;
  m_poly_force(i, 1) += f;
  m_poly_force(i, 2) += f;
}

void Polymer::add_force(int i, double fx, double fy, double fz)
{
  m_poly_force(i, 0) += fx;
  m_poly_force(i, 1) += fy;
  m_poly_force(i, 2) += fz;
}


// Input/Output function for read and write polymers
bool print_xyz(Polymer poly, std::string out_xyz )
{	
  // open stream to write xyz file
  std::ofstream output;
  output.open( out_xyz, std::fstream::app );

	// return error if read file fail
	if( output.fail() ) {
		throw "ERROR: Impossible to write xyz file "+out_xyz;
		return 1;
	}
	
  output << poly.get_poly_sphere() << std::endl << std::endl;
	for( int i = 0; i< poly.get_poly_sphere(); i++){
		output << "\tAu\t\t\t"<< poly.get_x(i) << "\t\t\t" << poly.get_y(i) << "\t\t\t" << poly.get_z(i) 
		 << "\t\t\t"<< poly.get_vx(i) << "\t\t\t" << poly.get_vy(i) << "\t\t\t" << poly.get_vz(i) << std::endl;
	}

  output.close();
  return 0;
}

arma::mat read_xyz(std::string in_xyz)
{
  double x, y, z,
         vx, vy, vz;
  arma::mat r_v(10, 6);
  
  std::string name;
  std::ifstream read_conf;
  read_conf.open( in_xyz, std::fstream::in );
  if( read_conf.fail() ){
		throw "ERROR: Impossible to open polymer coordinates file "+in_xyz;
	}
  
  std::string dummyLine;
  getline(read_conf, dummyLine);
  getline(read_conf, dummyLine);
  
  std::string line;
  int i=0;
   while(std::getline(read_conf, line)){
    std::stringstream iss(line);
    iss >> name >> x >> y >> z >> vx >> vy >> vz;
    r_v(i,0)=x;
    r_v(i,1)=y;
    r_v(i,2)=z;
    r_v(i,3)=vx;
    r_v(i,4)=vy;
    r_v(i,5)=vz;
    i++;
   }
  read_conf.close();
  return r_v;
}
