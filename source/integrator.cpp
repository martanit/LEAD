#include "integrator.h"

//integrator

void Integrator::langevin_overdamped()
{
  std::random_device rd;
  std::mt19937 mt(rd());
  std::normal_distribution<double> gauss_term(0., 1.);
  double stoch_term = std::sqrt(2.*m_integrator_gamma*m_integrator_temp/m_integrator_timestep);

  for(unsigned int i=0; i<(*m_poly).get_poly_sphere(); ++i){
    (*m_poly_new).set_x(pbc((*m_poly).get_x(i)+
                     ((*m_poly).get_fx(i)+
                     stoch_term*gauss_term(mt))*m_integrator_timestep/m_integrator_gamma), i);

    (*m_poly_new).set_y(pbc((*m_poly).get_y(i)+
                     ((*m_poly).get_fy(i)+
                     stoch_term*gauss_term(mt))*m_integrator_timestep/m_integrator_gamma), i);
   
    (*m_poly_new).set_z(pbc((*m_poly).get_z(i)+
                     ((*m_poly).get_fz(i)+
                     stoch_term*gauss_term(mt))*m_integrator_timestep/m_integrator_gamma), i);

    }
  m_poly=std::make_unique<Polymer>(*m_poly_new);
}

void Integrator::markov_chain()
{
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> transition_prob(0, 1);
    for(auto&& i: m_vector_extr){
       if((*i).get_r()!=((*m_poly).get_poly_sphere()-1))
             if ((( (*i).get_rate_r((int)(*i).get_r())) > transition_prob(mt)) 
                     and (*i).get_ctcf()[(*i).get_r()]!=1
                     and !(m_vector_extr.overlap_r(*i))) (*i).set_r((*i).get_r()+1); 
       if((*i).get_l()!=0)
             if ((( (*i).get_rate_l((*i).get_l())) > transition_prob(mt))  
                     and (*i).get_ctcf()[(*i).get_l()]!=(-1)
                     and !(m_vector_extr.overlap_l(*i))) (*i).set_l((*i).get_l()-1);
    }
}
