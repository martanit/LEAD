#include "polymer.h"

bool print_xyz(Polymer &, std::string);
bool read_xyz(Polymer &, std::string);

Polymer::Polymer(Parameters parm)
    : m_poly_nmonomers(parm.get_nmonomers()),
      m_poly_d(parm.get_diameter()), m_poly_spring(parm.get_spring()),
      uniform01(0., 1.){
  this->set_size();
  this->reset_force();
  this->poly_configuration();
}

Polymer::Polymer(Parameters parm, std::string poly_xyz)
    : m_poly_nmonomers(parm.get_nmonomers()),
      m_poly_d(parm.get_diameter()), m_poly_spring(parm.get_spring()),
      uniform01(0., 1.) {
  this->set_size();
  this->reset_force();
  read_xyz((*this), poly_xyz);
  // this->control_poly;
}

Polymer::~Polymer() {}

void Polymer::set_size() {
  m_poly_x.resize(m_poly_nmonomers);
  m_poly_y.resize(m_poly_nmonomers);
  m_poly_z.resize(m_poly_nmonomers);
  m_poly_fx.resize(m_poly_nmonomers);
  m_poly_fy.resize(m_poly_nmonomers);
  m_poly_fz.resize(m_poly_nmonomers);
}

void Polymer::first_monomer() {
  m_poly_x[0] = (uniform01(mt)); // x
  m_poly_y[0] = (uniform01(mt)); // y
  m_poly_z[0] = (uniform01(mt)); // z
}

void Polymer::poly_configuration() {
  unsigned int stuck;
  double d = m_poly_d;
  this->first_monomer();
  for (unsigned int i = 1; i < m_poly_nmonomers; i++) {
    stuck = 0;
    bool is_overlap = true;
    while (is_overlap) {
      stuck++;
      double t = 2 * M_PI * uniform01(mt);
      // uniform distribution on monomers
      double p = acos(1 - 2 * uniform01(mt));

      m_poly_x[i] = (m_poly_x[i - 1] + d * std::sin(p) * std::cos(t));
      m_poly_y[i] = (m_poly_y[i - 1] + d * std::sin(p) * std::sin(t));
      m_poly_z[i] = (m_poly_z[i - 1] + d * std::cos(p));
      is_overlap = this->is_overlap(i);
	if(stuck > 10000){
	       std::cout << "I'm not able to construct polymer (due to randomness). Please retry!" << std::endl;
       	       exit(1);
	}
    }
  }
}

bool Polymer::is_overlap(int lenght) {
  bool is_overlap = false;
  for (unsigned int j = 0; j < lenght; j++) {
    if (this->dist(j, lenght) < m_poly_d) {
      is_overlap = true;
      break;
    } else
      is_overlap = false;
  }
  return is_overlap;
}

double Polymer::dist(const int &i, const int &j) {

  d = std::sqrt(std::pow(m_poly_x[i] - m_poly_x[j], 2) +
                std::pow(m_poly_y[i] - m_poly_y[j], 2) +
                std::pow(m_poly_z[i] - m_poly_z[j], 2));
  return (d < 10E-9) ? 10E-9 : d;
}

void Polymer::center() {
	this->set_cm();
	std::for_each(m_poly_x.begin(), m_poly_x.end(), [this]
			(auto& x) { x -= x_cm; });
	std::for_each(m_poly_y.begin(), m_poly_y.end(), [this]
			(auto& y) { y -= y_cm; });
	std::for_each(m_poly_z.begin(), m_poly_z.end(), [this]
			(auto& z) { z -= z_cm; });
}

void Polymer::set_cm() {
	sum_x = 0;
	sum_y = 0;
	sum_z = 0;

	for(unsigned int i=0; i<m_poly_nmonomers; ++i){
		sum_x += get_x(i);
		sum_y += get_y(i);
		sum_z += get_z(i);
	}
	x_cm = sum_x/m_poly_nmonomers;
	y_cm = sum_y/m_poly_nmonomers;
	z_cm = sum_z/m_poly_nmonomers;
}

void Polymer::reset_force() {
  std::fill(m_poly_fx.begin(), m_poly_fx.end(), 0.);
  std::fill(m_poly_fy.begin(), m_poly_fy.end(), 0.);
  std::fill(m_poly_fz.begin(), m_poly_fz.end(), 0.);
}

void Polymer::add_force(const int &i, const double &fx, const double &fy,
                        const double &fz) {
  m_poly_fx[i] += fx;
  m_poly_fy[i] += fy;
  m_poly_fz[i] += fz;
}

void Polymer::reset_energy() { m_poly_e = 0.; }

void Polymer::add_energy(const double &e) { m_poly_e += e; }

// Input/Output function for read and write polymers
bool print_xyz(Polymer &poly, std::string out_xyz) {
  // open stream to write xyz file
  std::ofstream output;
  output.open(out_xyz, std::fstream::app);

  // return error if read file fail
  if (output.fail()) {
    throw "ERROR: Impossible to write xyz file " + out_xyz;
    return 1;
  }

  output << poly.m_poly_nmonomers << std::endl << std::endl;
  for (int i = 0; i < poly.m_poly_nmonomers; i++) {
    output << "\tAu" << "\t\t\t" 
           << poly.m_poly_x[i] << "\t\t\t"
           << poly.m_poly_y[i] << "\t\t\t" 
           << poly.m_poly_z[i] << std::endl;
  }

  output.close();
  return 0;
}

bool read_xyz(Polymer &poly, std::string in_xyz) {
  std::string name;
  std::ifstream read_conf;
  read_conf.open(in_xyz, std::fstream::in);
  if (read_conf.fail()) {
    throw "ERROR: Impossible to open polymer coordinates file " + in_xyz;
    return 1;
  }

  std::string dummyLine;
  getline(read_conf, dummyLine);
  getline(read_conf, dummyLine);

  std::string line;
  int i = 0;
  while (std::getline(read_conf, line)) {
    std::stringstream iss(line);
    iss >> name >> poly.m_poly_x[i] >> poly.m_poly_y[i] >> poly.m_poly_z[i];
    i++;
  }
  read_conf.close();
  return 0;
}
