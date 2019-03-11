#ifndef POLYMER_H_
#define POLYMER_H_

#include "parameters.h"
#include <cmath>
#include <vector>
#include <random>
#include <fstream>

class Polymer {
public:
	// default constructor: assing default parameters to polymer
	Polymer();
	// assign passed parameters to polymer
	Polymer(std::string);
	// contruct polymer from x, y, z coordinate file and assign parameters
	Polymer(std::vector<double>, std::vector<double>, std::vector<double>, std::string);
	~Polymer();
	
	
	// place first sphere of polymer
	void first_sphere();
	// place randomly other sphere
	void poly_configuration();
	// check if spheres overlap (self avoiding polymer)
	bool is_overlap(int );
	// print coordinates of polymer in vmd like .xyz coordinates file
	bool print_xyz(std::string );

	// auxiliary function
	int get_nsphere();

private:
	std::vector<double> m_polX, m_polY, m_polZ; //vector of x,y,z coordinates of polymer
	int m_psphere;
	Parameters m_parm;
};

#endif /* POLYMER_H_ */
