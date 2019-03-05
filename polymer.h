/*
 * Polymer.h
 *
 *  Created on: Feb 19, 2019
 *      Author: martanit
 */

#ifndef POLYMER_H_
#define POLYMER_H_

#include "Parameters.h"
#include <cmath>

class Polymer {
public:
	Polymer();
	Polymer(
	Polymer(std::array<double, int >, std::array<double, int >, std::array<double, int >);
	~Polymer();
	
	

	void first_sphere();
	void poly_configuration();
	bool is_overlap(int &);

	int get_nsphere();
	
private:
	std::array<double, int > m_polX, m_polY, m_polZ; //array of x,y,z coordinates of polymer
	Parameters m_parm;

};

#endif /* POLYMER_H_ */
