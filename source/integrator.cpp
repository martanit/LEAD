#include "integrator.h"

//integrator

void Integrator::langevin_overdamped()
{
  for(unsigned int i=0; i<(*m_poly).get_poly_sphere(); ++i){
    (*m_poly).set_x(pbc((*m_poly).get_x(i)+
                     ((*m_poly).get_fx(i)+
                     stoch_term*gauss_term(mt))*m_integrator_timestep/m_integrator_gamma), i);

    (*m_poly).set_y(pbc((*m_poly).get_y(i)+
                     ((*m_poly).get_fy(i)+
                     stoch_term*gauss_term(mt))*m_integrator_timestep/m_integrator_gamma), i);
   
    (*m_poly).set_z(pbc((*m_poly).get_z(i)+
                     ((*m_poly).get_fz(i)+
                     stoch_term*gauss_term(mt))*m_integrator_timestep/m_integrator_gamma), i);
    
    (*m_poly).set_vx((1-(m_integrator_gamma/(*m_poly).get_poly_mass())*
                      m_integrator_timestep)*(*m_poly).get_vx(i)+
                      std::sqrt(2*m_integrator_temp*m_integrator_timestep*
                     m_integrator_gamma/(*m_poly).get_poly_mass())*
                     gauss_term(mt), i);

    (*m_poly).set_vy((1-(m_integrator_gamma/(*m_poly).get_poly_mass())*
                      m_integrator_timestep)*(*m_poly).get_vy(i)+
                      std::sqrt(2*m_integrator_temp*m_integrator_timestep*
                     m_integrator_gamma/(*m_poly).get_poly_mass())*
                     gauss_term(mt), i);
    
    (*m_poly).set_vz((1-(m_integrator_gamma/(*m_poly).get_poly_mass())*
                      m_integrator_timestep)*(*m_poly).get_vz(i)+
                      std::sqrt(2*m_integrator_temp*m_integrator_timestep*
                     m_integrator_gamma/(*m_poly).get_poly_mass())*
                     gauss_term(mt), i);

  }
}

void Integrator::markov_chain()
{
    for(auto&& i: m_vector_extr){
        // right side go foward to right
       if((*i).get_r()!=((*m_poly).get_poly_sphere()-1))
             if (( (*i).get_rate_fwr() > transition_prob(mt)) 
                     and ((*i).get_ctcf()[(*i).get_r()]!=1 or permeability_ctcf < transition_prob(mt)) 
                     and !(m_vector_extr.overlap_rl(*i))) (*i).set_r((*i).get_r()+1); 
       // right side go backward to left
             if (( (*i).get_rate_bwr() > transition_prob(mt)) 
                     and ((*i).get_ctcf()[(*i).get_r()]!=(-1) or permeability_ctcf < transition_prob(mt)) 
                     and !(m_vector_extr.overlap_rr(*i))) (*i).set_r((*i).get_r()-1); 
       
       // left side go forward to left
       if((*i).get_l()!=0)
             if (( (*i).get_rate_fwl() > transition_prob(mt))  
                     and ((*i).get_ctcf()[(*i).get_l()]!=(-1) or permeability_ctcf < transition_prob(mt)) 
                     and !(m_vector_extr.overlap_lr(*i))) (*i).set_l((*i).get_l()-1);
       // left side go backward to right
             if (( (*i).get_rate_bwl() > transition_prob(mt))  
                     and ((*i).get_ctcf()[(*i).get_l()]!=1 or permeability_ctcf < transition_prob(mt)) 
                     and !(m_vector_extr.overlap_ll(*i))) (*i).set_l((*i).get_l()+1);
    }
}

void Integrator::markov_chain_rate()
{
    for(auto&& i: m_vector_extr){
       if((*i).get_r()!=((*m_poly).get_poly_sphere()-1))
             if (( (*i).get_rate_vfwr((int)(*i).get_r()) > transition_prob(mt)) 
                     and (*i).get_ctcf()[(*i).get_r()]!=1
                     and !(m_vector_extr.overlap_rl(*i))) (*i).set_r((*i).get_r()+1); 

       if (( (*i).get_rate_vbwr((int)(*i).get_r()) > transition_prob(mt)) 
                     and (*i).get_ctcf()[(*i).get_r()]!=(-1)
                     and !(m_vector_extr.overlap_rr(*i))) (*i).set_r((*i).get_r()-1); 
       
       if((*i).get_l()!=0)
             if (( (*i).get_rate_vfwl((int)(*i).get_l()) > transition_prob(mt))  
                     and (*i).get_ctcf()[(*i).get_l()]!=(-1)
                     and !(m_vector_extr.overlap_lr(*i))) (*i).set_l((*i).get_l()-1);
       
       if (( (*i).get_rate_vbwl((int)(*i).get_l()) > transition_prob(mt))  
                     and (*i).get_ctcf()[(*i).get_l()]!=1
                     and !(m_vector_extr.overlap_ll(*i))) (*i).set_l((*i).get_l()+1);
    }
}
