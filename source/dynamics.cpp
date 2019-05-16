#include "dynamics.h"

void Dynamics::run()
{
  for(unsigned int i=0; i<m_dynamics_nstep; ++i){ 
      
      (*m_poly).reset_force();
      
      Potential::set_new_polymer(*m_poly);
      this->lennard_jones_f();
      this->harmonic_spring_f();
      m_poly = std::make_unique<Polymer>(Potential::get_poly());
      
      Integrator::set_new_polymer(*m_poly);
      this->euler();
      
      m_poly = std::make_unique<Polymer>(Integrator::get_poly());  
      if(i%m_dynamics_print == 0) print_xyz(*m_poly, "traj.xyz");
    
  }
}
