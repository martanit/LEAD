/*
 * main.cpp
 *
 *  Created on: Feb 26, 2019
 *      Author: martanit
 */

#include "dynamics.h"
#include "extruder.h"
#include "integrator.h"
#include "parameters.h"
#include "polymer.h"
#include "potential.h"
#include "vector_extruder.h"

int main() {
  Parameters parm("input/parameters.in", "input/ctcf.in",
                  "input/coupling_probability.in");
  Polymer poly_init(parm);
  Extruder extr(parm);
  VectorExtruder v_extr(parm, extr, poly_init);

  Dynamics dyn(poly_init, v_extr, parm);

  print_xyz(poly_init, "output/traj.xyz");

  dyn.run_extrusion();

  Polymer poly_last = dyn.get_poly();
  print_xyz(poly_last, "output/traj.xyz");

  return 0;
}
