/*
 * main.cpp
 *
 *  Created on: Feb 26, 2019
 *      Author: martanit
 */

#include "parameters.h"
#include "polymer.h"
#include "potential.h"
#include "dynamics.h"

int main() {
  
  Parameters par("parameters.in");
  Polymer poly_init(par, "initial_chain.xyz");
  Dynamics dyn(poly_init, par);
  dyn.run(); 
	return 0;
}
