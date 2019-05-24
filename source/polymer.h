#ifndef POLYMER_H_
#define POLYMER_H_

#include "parameters.h"
#include "utils.h"
#include <cmath>
#include <random>
#include <fstream>
#include <vector>

class Polymer {
  
  public:
  	// default constructor: assing default parameters to polymer
  	Polymer();
  	// assign passed parameters to polymer
  	Polymer(Parameters );
  	// contruct polymer from x, y, z coordinate file and assign parameters
  	Polymer(Parameters, std::string);
  	
    ~Polymer();
    	
	void set_size();
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
    inline const double& get_x(int i) const { return m_poly_x[i]; };
    inline const double& get_y(int i) const { return m_poly_y[i]; };
    inline const double& get_z(int i) const { return m_poly_z[i]; };
    
    inline const double& get_vx(int i) const { return m_poly_vx[i]; };
    inline const double& get_vy(int i) const { return m_poly_vy[i]; };
    inline const double& get_vz(int i) const { return m_poly_vz[i]; };
   
    inline const double& get_fx(int i) const { return m_poly_fx[i]; };
    inline const double& get_fy(int i) const { return m_poly_fy[i]; };
    inline const double& get_fz(int i) const { return m_poly_fz[i]; };
    
    inline void set_x(double x, int i) { m_poly_x[i] = x; };
    inline void set_y(double y, int i) { m_poly_y[i] = y; };
    inline void set_z(double z, int i) { m_poly_z[i] = z; };
    
    inline void set_vx(double vx, int i) { m_poly_vx[i] = vx; };
    inline void set_vy(double vy, int i) { m_poly_vy[i] = vy; };
    inline void set_vz(double vz, int i) { m_poly_vz[i] = vz; };
    
    void reset_force();    
    void set_force(int, double, double, double);
    void add_force(int, double, double, double);    
    
    // function to access polymer parameters
    inline const int& get_poly_sphere() const { return m_poly_sphere; };
    inline const float& get_poly_dist() const { return m_poly_dist; };
    inline const float& get_bond() const  { return m_poly_bond; };
    inline const float& get_poly_mass() const { return m_poly_mass; };

    friend bool print_xyz(Polymer&, std::string);
    friend bool read_xyz(Polymer&, std::string);
  
  private:
    
    // polymer parameters
	float m_poly_mass = 1.;		    // mass
	int m_poly_sphere = 100;		// number of interacting sphere
	float m_poly_dist = 2.;		    // length of bond [angstrom]
	float m_poly_bond = 10.;		// spring constant of harmonic oscillator (bond)
	float m_poly_hradius = 0.2;		// hard core radius of sphere
    std::string m_poly_init = "init.dat";
    
    std::vector<double> m_poly_x; 
    std::vector<double> m_poly_y; 
    std::vector<double> m_poly_z; 
    
    std::vector<double> m_poly_vx; 
    std::vector<double> m_poly_vy; 
    std::vector<double> m_poly_vz; 
    
    std::vector<double> m_poly_fx; 
    std::vector<double> m_poly_fy; 
    std::vector<double> m_poly_fz; 
    
    Utils m_conf;  
};


#endif /* POLYMER_H_ */
