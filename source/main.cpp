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
#include "extruder.h"

int main() {
  
  Parameters parm("parameters.in", "ctcf.in", "coupling_probability.in");
    
  //Polymer poly_init(poly_par, "initial_chain.xyz");
  Polymer poly_init(parm);
  Extruder extr(parm, poly_init);
  Dynamics dyn(poly_init, extr, parm);
  print_xyz(poly_init, "traj.xyz");
  
  dyn.run(); 
  
  Polymer poly_last = dyn.get_poly();
  print_xyz(poly_last, "traj.xyz");
  
  return 0;
}
