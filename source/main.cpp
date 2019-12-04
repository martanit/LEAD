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
#include "cohesin_field.h"

#include<chrono>
#include<iostream>
#include<cstring>

int main(int argc, char** argv) {

  bool compute_energy=false;
  bool extrusion=false;
  bool rouse=false;
  bool soft_core=false;

  std::string parm_input;
  std::string parm_output;
  std::string traj_output;

  if(argc<7){
	  std::cerr << "Usage: -i parm_input -o parm_output -trj trj_output" <<std::endl;
	  return 1;
  }
  for(int idx=1; idx<6; idx=idx+2){
	  if(std::string(argv[idx]) == "-i"){
		  parm_input = std::string(argv[idx+1]);
		  std::cout << "Input parameters read from: " << parm_input << std::endl;
	  }
	  else if(std::string(argv[idx]) == "-o"){
	  	parm_output = std::string(argv[idx+1]);
		std::cout << "Parameters used written on: " << parm_output+".out" << std::endl;
	  }
	  else if(std::string(argv[idx]) == "-trj"){
  		traj_output = std::string(argv[idx+1]);
		std::cout << "Trajectory written on: " << traj_output+".xyz" << std::endl;
	  }
 	  else { 
		std::cerr << "Usage: -i parm_input -o parm_output -trj trj_output" <<std::endl;
	  	return 1;
	  }
  }
  if(argc == 7) std::cout << "Using Lennard Jones potential for spring chain without loop extrusion" << std::endl;
  if(argc > 7){
  for(int idx=7; idx<argc; idx++){ 
	  if(std::string(argv[idx])== "-energy"){
	  	compute_energy=true;
		std::cout <<"Computing energy..." << std::endl;
	  }
	  else if(std::string(argv[idx])== "-le"){
	  	extrusion=true;
		std::cout << "Loop extrusion activated" << std::endl;
	  }
  	  else if(std::string(argv[idx])== "-rouse"){
	  	rouse=true;
		std::cout << "Using only Rouse polymer" << std::endl;
	  }
  	  else if(std::string(argv[idx])== "-soft_core"){
	  	  soft_core=true;
		  std::cout << "Using only Rouse polymer with soft core repulsion" << std::endl;  
	  }
	  else{
		  std::cerr << "Unknown flag, options are: -energy, -le, -rouse, -soft_core" << std::endl;
		  return 1;
  	}
  	}
  }
  
auto begin = std::chrono::high_resolution_clock::now();

  if(!extrusion){
  Parameters parm(parm_input, parm_output+".out");
  Polymer poly_init(parm);
  poly_init.center();
  print_xyz(poly_init, traj_output+".xyz");
  
  Dynamics dyn(poly_init, parm);
 
  dyn.run(rouse, soft_core, compute_energy, traj_output+".xyz");
 
  Polymer poly_last = dyn.get_poly();
  print_xyz(poly_last, traj_output+".xyz");
  
  }
  else{
  
  Parameters parm(parm_input, parm_output+".out", "input/ctcf.in",
                  "input/coupling_probability.in");  
  Polymer poly_init(parm);
  poly_init.center();
  print_xyz(poly_init, traj_output+".xyz");
  
  Extruder extr(parm);
  VectorExtruder v_extr(parm, extr, poly_init);
  Dynamics dyn(poly_init, v_extr, parm);
 
  dyn.run_extrusion(rouse, soft_core, compute_energy, traj_output+".xyz");
  
  Polymer poly_last = dyn.get_poly();
  print_xyz(poly_last, traj_output+".xyz");
  
  }
 auto end = std::chrono::high_resolution_clock::now();
  std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1e6 << "ms\n";


  return 0;
}
