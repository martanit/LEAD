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
  
  Polymer_Parameters poly_par("parameters.in");
  Potential_Parameters pot_par("parameters.in");
  Dynamics_Parameters dyn_par("parameters.in");
  
  //Polymer poly_init(poly_par, "initial_chain.xyz");
  Polymer poly_init(poly_par);
  Potential pot(poly_init, pot_par); 
  Dynamics dyn(poly_init, dyn_par, pot);

  print_xyz(poly_init, "traj.xyz");
  dyn.run(); 
  
  Polymer poly_last = dyn.get_poly();
  
  print_xyz(poly_last, "traj.xyz");
  
  return 0;
}
