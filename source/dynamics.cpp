#include "dynamics.h"

void Dynamics::run_extrusion()
{  
    bool compute_energy = false;
    for(unsigned int i=0; i<m_dynamics_nstep; ++i){ 
        if(i%m_dynamics_print == 0)    
            m_vector_extr.update(*m_poly);
      
        (*m_poly).reset_force();
       
        Potential::set_new_polymer(*m_poly);
        Potential::set_new_extruder(m_vector_extr);

        this->extruder_spring_f();
        this->lennard_jones_f(i, true, compute_energy);
        this->harmonic_spring_f();
        
        m_poly = std::make_unique<Polymer>(Potential::get_poly());
        m_vector_extr = Potential::get_extr();

        Integrator::set_new_polymer(*m_poly);
        Integrator::set_new_extruder(m_vector_extr);
      
        this->markov_chain();
        this->langevin_overdamped();

        m_poly = std::make_unique<Polymer>(Integrator::get_poly()); 
        m_vector_extr = Integrator::get_extr();
        
        if(i%m_dynamics_print == 0) {
            print_xyz(*m_poly, "output/traj.xyz");
            int num_extr=0;
            for (auto &i : m_vector_extr){
              print_r(*m_poly, *i, "output/loop_extrusion_"+std::to_string(num_extr)+".r");
              ++num_extr;
            }
        }
    }
}

