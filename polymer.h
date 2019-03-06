#ifndef POLYMER_H_
#define POLYMER_H_

#include "parameters.h"
#include <cmath>
#include <vector>
#include <random>

class Polymer {
public:
	Polymer();
	Polymer(std::vector<double>, std::vector<double>, std::vector<double>);
	~Polymer();
	
	

	void first_sphere();
	void poly_configuration();
	bool is_overlap(int );
	
	int get_nsphere();
	
private:
	std::vector<double> m_polX, m_polY, m_polZ; //vector of x,y,z coordinates of polymer
	int m_psphere;
	Parameters m_parm=Parameters("data.dat");
};

#endif /* POLYMER_H_ */
