#ifndef POLYMER_H_
#define POLYMER_H_

#include "parameters.h"
#include <cmath>
#include <vector>
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

    // print and read coordinates of polymer in vmd like .xyz coordinates file
	  bool print_xyz(std::string );
    bool read_xyz(std::string);

    // calculate distance between atoms
    double dist(int , int );
    // set periodic boundary conditions
    double pbc(double );  
    
    // fuction to access sphere coordinates
    double get_x(int i){ return m_poly_r_v(i,0); };
    double get_y(int i){ return m_poly_r_v(i,1); };
    double get_z(int i){ return m_poly_r_v(i,2); };
  
  private:
    arma::mat m_poly_r_v; //vector of x, y, z, vx, vy, vz
	
    // polymer parameters
	  float m_poly_mass = 1.;		// mass
	  int m_poly_sphere = 100;		// number of interacting sphere
	  float m_poly_dist = 10.;		// length of bond [angstrom]
	  float m_poly_bond = 15.;		// spring constant of harmonic oscillator (bond)
	  float m_poly_hradius = 5.;		// hard core radius of sphere
};

#endif /* POLYMER_H_ */
