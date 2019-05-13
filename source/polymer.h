#ifndef POLYMER_H_
#define POLYMER_H_

#include "parameters.h"
#include "utils.h"
#include <cmath>
#include <random>
#include <fstream>
#include <armadillo>

class Polymer {
  
  public:
  	// default constructor: assing default parameters to polymer
  	Polymer();
  	// assign passed parameters to polymer
  	Polymer(Parameters );
  	// contruct polymer from x, y, z coordinate file and assign parameters
  	Polymer(Parameters, std::string);
  	~Polymer();
    	
	  // place first sphere of polymer
	  void first_sphere();
	  // place randomly other sphere
	  void poly_configuration();
	  // check if spheres overlap (self avoiding polymer)
	  bool is_overlap(int );
	  // place randomly velocities with center of mass
    // centered in zero
    void poly_velocity();
    
    // calculate distance between atoms
    double dist(int , int );

    // fuction to access sphere coordinates, velocities and forces
    const double& get_x(int i) const { return m_poly_r_v(i,0); };
    const double& get_y(int i) const { return m_poly_r_v(i,1); };
    const double& get_z(int i) const { return m_poly_r_v(i,2); };
    
    const double& get_vx(int i) const { return m_poly_r_v(i,3); };
    const double& get_vy(int i) const { return m_poly_r_v(i,4); };
    const double& get_vz(int i) const { return m_poly_r_v(i,5); };
    
    const double& get_force_x(int i) const {return m_poly_force(i,0); }; 
    const double& get_force_y(int i) const {return m_poly_force(i,1); }; 
    const double& get_force_z(int i) const {return m_poly_force(i,2); }; 

    void set_x(double x, int i) { m_poly_r_v(i,0) = x; };
    void set_y(double y, int i) { m_poly_r_v(i,1) = y; };
    void set_z(double z, int i) { m_poly_r_v(i,2) = z; };
    
    void set_vx(double vx, int i) { m_poly_r_v(i,3) = vx; };
    void set_vy(double vy, int i) { m_poly_r_v(i,4) = vy; };
    void set_vz(double vz, int i) { m_poly_r_v(i,5) = vz; };
    void reset_force();    
    void set_force(int, double, double, double);
    void add_force(int, double, double, double);    
    // overload of add force: for spring chain only f=fx=fy=fz needed
    void add_force(int, double);
    
    // function to access polymer parameters
    const int& get_poly_sphere() const { return m_poly_sphere; }
    const float& get_poly_dist() const { return m_poly_dist; }
    const float& get_bond() const  { return m_poly_bond; }
    const float& get_poly_mass() const { return m_poly_mass; }

  private:
    arma::mat m_poly_r_v; //vector of x, y, z, vx, vy, vz
    arma::mat m_poly_force;  
    Utils m_conf;  
    
    // polymer parameters
	  float m_poly_mass = 1.;		// mass
	  int m_poly_sphere = 100;		// number of interacting sphere
	  float m_poly_dist = 2.;		// length of bond [angstrom]
	  float m_poly_bond = 10.;		// spring constant of harmonic oscillator (bond)
	  float m_poly_hradius = 0.2;		// hard core radius of sphere
};

bool print_xyz(Polymer, std::string );
arma::mat read_xyz(std::string);

#endif /* POLYMER_H_ */
