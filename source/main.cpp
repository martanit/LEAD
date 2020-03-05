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

#include<chrono>
#include<iostream>
#include<cstring>
    
bool print_box(Polymer, Parameters, std::string);
int main(int argc, char** argv) {

    bool compute_energy=false;
    bool expl_extrusion=false;
    bool eff_extrusion=false;
    bool rouse=false;
    bool soft_core=false;
    bool lennard_jones=false;

    std::string parm_input;
    std::string parm_output;
    std::string traj_output;

    if(argc<8) {
        std::cerr << "Usage: -i parm_input -o parm_output -trj trj_output --POTENTIAL" <<std::endl;
        std::cerr << "Choices for POTENTIAL are --rouse, --soft-core, --lennard-jones" <<std::endl;
        return 1;
    }
    for(int idx=1; idx<6; idx=idx+2) {
        if(std::string(argv[idx]) == "-i") {
            parm_input = std::string(argv[idx+1]);
            std::cerr << "Input parameters read from: " << parm_input << std::endl;
        }
        else if(std::string(argv[idx]) == "-o") {
            parm_output = std::string(argv[idx+1]);
            std::cerr << "Parameters used written on: " << parm_output << ".out" << std::endl;
        }
        else if(std::string(argv[idx]) == "-trj") {
            traj_output = std::string(argv[idx+1]);
            std::cerr << "Trajectory written on: " << traj_output << ".xyz" << std::endl;
        }
        else {
            std::cerr << "Usage: -i parm_input -o parm_output -trj trj_output" <<std::endl;
            return 1;
        }
    }
    for(int idx=7; idx<argc; idx++) {
        if(std::string(argv[idx])== "--energy") {
            compute_energy=true;
            std::cerr <<"Computing energy..." << std::endl;
        }
        else if(std::string(argv[idx])== "--le-expl") {
            expl_extrusion=true;
            std::cerr << "Loop extrusion with explicit cohesin activated" << std::endl;
        }
        else if(std::string(argv[idx])== "--le-eff") {
            eff_extrusion=true;
            std::cerr << "Loop extrusion with effective potential activated" << std::endl;
        }
        else if(std::string(argv[idx])== "--rouse") {
            rouse=true;
            std::cerr << "Using Rouse polymer" << std::endl;
        }
        else if(std::string(argv[idx])== "--soft-core") {
            soft_core=true;
            std::cerr << "Using Rouse polymer with soft core repulsion" << std::endl;
        }
        else if(std::string(argv[idx])== "--lennard-jones") {
            lennard_jones=true;
            std::cerr << "Using Rouse polymer with Lennard Jones interaction" << std::endl;
        }
        else {
            std::cerr << "Unknown flag, options are: --energy, --le-expl, le-eff, --rouse, --soft-core, --lennard-jones" << std::endl;
            return 1;
        }
    }
    if((rouse==true and soft_core==true) or (rouse==true and lennard_jones==true) or (soft_core==true and lennard_jones==true)) {
        std::cerr << "ERROR: you have to use only one potential" << std::endl;
        return 1;
    }
    if(expl_extrusion==true and eff_extrusion==true) {
        std::cerr << "ERROR: you have to choose the explicit OR the effective model" << std::endl;
        return 1;
    }

    auto begin = std::chrono::high_resolution_clock::now();

    try{    
    if(!eff_extrusion and !expl_extrusion) {
        Parameters parm(parm_input, parm_output+".out", "no_extr");
        Polymer poly_init(parm);
        print_box(poly_init,parm, parm_output+".box");
        Dynamics dyn(poly_init, parm);

        dyn.run(rouse, soft_core, lennard_jones, compute_energy, traj_output+".xyz");
    }
    else if(eff_extrusion) {
        Parameters parm(parm_input, parm_output+".out", "eff_extrusion");
        Polymer poly_init(parm);
        print_box(poly_init,parm, parm_output+".box");
        Dynamics dyn(poly_init, parm);

        dyn.run_effective(rouse, soft_core, lennard_jones, compute_energy, traj_output+".xyz");
    }
    else if(expl_extrusion){
        Parameters parm(parm_input, parm_output+".out", "expl_extrusion");
        Polymer poly_init(parm);
        print_box(poly_init, parm, parm_output+".box");
        Extruder extr(parm);
        VectorExtruder v_extr(parm, extr, poly_init);
        Dynamics dyn(poly_init, v_extr, parm);

        dyn.run_explicit(rouse, soft_core, lennard_jones, compute_energy, traj_output+".xyz");
    }
    }
    catch(std::string s){
 	 std::cerr <<  s << '\n';	
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::cerr << std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1e6 << "ms\n";

    return 0;
}

bool print_box(Polymer poly, Parameters parm, std::string string) {
    
    std::ofstream set_box;
    set_box.open(string, std::ofstream::out);
    // return error if read file fail
    if (set_box.fail()) {
        throw "ERROR: Impossible to write box to file "+ string;
        return 1;
    }

    set_box << poly.get_xcm()-parm.get_box_length()/2. << " " << poly.get_xcm()+parm.get_box_length()/2. << std::endl;
    set_box << poly.get_ycm()-parm.get_box_length()/2. << " " << poly.get_ycm()+parm.get_box_length()/2. << std::endl;
    set_box << poly.get_zcm()-parm.get_box_length()/2. << " " << poly.get_zcm()+parm.get_box_length()/2. << std::endl;
    
    return 0;
}
