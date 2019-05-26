#include "dynamics.h"

void Dynamics::run()
{
    for(unsigned int i=0; i<m_dynamics_nstep; ++i){ 
      
        if(i%m_dynamics_nstep == 0){
            
        
        
        
        
        } 
        
        (*m_poly).reset_force();
      
        Potential::set_new_polymer(*m_poly);
        Potential::set_new_extruder(m_extr);
      
        this->extruder_spring_f();
        this->lennard_jones_f();
        this->harmonic_spring_f();
        
        m_poly = std::make_unique<Polymer>(Potential::get_poly());
	for (auto &i : Potential::get_extr())
		m_extr.push_back(std::make_unique<Extruder>(i)); 

        Integrator::set_new_polymer(*m_poly);
        Integrator::set_new_extruder(m_extr);
      
        this->markov_chain();
        this->langevin_overdamped();
        
        m_poly = std::make_unique<Polymer>(Integrator::get_poly()); 
	for (const auto &i : Integrator::get_extr())
		m_extr.push_back(std::make_unique<Extruder>(i)); 

        if(i%m_dynamics_print == 0) {
            print_xyz(*m_poly, "output/traj.xyz");
     //       print_r(*m_poly, *m_extr, "output/loop_extrusion.r");
        }
    }
}
