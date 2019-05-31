#include "integrator.h"

//integrator

void Integrator::langevin_overdamped()
{
//#pragma omp parallel for  
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
  }
}

void Integrator::markov_chain()
{
    for(auto&& i: m_vector_extr){
       if((*i).get_r()!=((*m_poly).get_poly_sphere()-1))
             if (( (*i).get_rate_r() > transition_prob(mt)) 
                     and ((*i).get_ctcf()[(*i).get_r()]!=1 or permeability_ctcf < transition_prob(mt)) 
                     and !(m_vector_extr.overlap_r(*i))) (*i).set_r((*i).get_r()+1); 
       if((*i).get_l()!=0)
             if (( (*i).get_rate_l() > transition_prob(mt))  
                     and ((*i).get_ctcf()[(*i).get_l()]!=(-1) or permeability_ctcf < transition_prob(mt)) 
                     and !(m_vector_extr.overlap_l(*i))) (*i).set_l((*i).get_l()-1);
    }
}

void Integrator::markov_chain_rate()
{
    for(auto&& i: m_vector_extr){
       if((*i).get_r()!=((*m_poly).get_poly_sphere()-1))
             if (( (*i).get_rate_vr((int)(*i).get_r()) > transition_prob(mt)) 
                     and (*i).get_ctcf()[(*i).get_r()]!=1
                     and !(m_vector_extr.overlap_r(*i))) (*i).set_r((*i).get_r()+1); 
       if((*i).get_l()!=0)
             if (( (*i).get_rate_vl((int)(*i).get_l()) > transition_prob(mt))  
                     and (*i).get_ctcf()[(*i).get_l()]!=(-1)
                     and !(m_vector_extr.overlap_l(*i))) (*i).set_l((*i).get_l()-1);
    }
}
