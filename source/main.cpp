/*
 * main.cpp
 *
 *  Created on: Feb 26, 2019
 *      Author: martanit
 */

#include "parameters.h"
#include "polymer.h"
#include "potential.h"

int main() {
	
	Polymer poly2("data.dat");
	poly2.poly_configuration();
	poly2.print_xyz( "prova.xyz");
  Potential pot(poly2);
  std::vector<double> vec;
  std::vector<double> vec2;
  vec = pot.lennard_jones(4);
  vec2 = pot.harmonic_spring();

  for(int i=0; i<vec.size(); i++){
     std::cout<< vec[i]<<std::endl;
  }
  std::cout<< "**" << std::endl;
  for(int i=0; i<vec2.size(); i++){
     std::cout<< vec2[i]<<std::endl;
  }
  
  
	return 0;
}
