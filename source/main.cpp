/*
 * main.cpp
 *
 *  Created on: Feb 26, 2019
 *      Author: martanit
 */

#include "parameters.h"
#include "polymer.h"
#include "potential.h"
#include "integrator.h"
#include "dynamics.h"

int main() {
  
  Parameters parm("parameters.in");
  
 Polymer poly_init(parm, "initial_chain_100.xyz");
//  Polymer poly_init(parm);
  Dynamics dyn(poly_init, parm);
  print_xyz(poly_init, "traj.xyz");
  
  dyn.run(); 
  
  Polymer poly_last = dyn.get_poly();
  print_xyz(poly_last, "traj.xyz");
  
  return 0;
}
