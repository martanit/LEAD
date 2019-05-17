#include "integrator.h"

//integrator

void Integrator::euler()
{
  for(unsigned int i=0; i<(*m_poly).get_poly_sphere(); ++i){
  
    (*m_poly_new).set_x((*m_poly).get_x(i)+(*m_poly).get_vx(i)*m_integrator_timestep, i);
    (*m_poly_new).set_y((*m_poly).get_y(i)+(*m_poly).get_vy(i)*m_integrator_timestep, i);
    (*m_poly_new).set_z((*m_poly).get_z(i)+(*m_poly).get_vz(i)*m_integrator_timestep, i);
    
    (*m_poly_new).set_vx((*m_poly).get_vx(i)+(*m_poly).get_fx(i)/(*m_poly).get_poly_mass()*m_integrator_timestep, i);
    (*m_poly_new).set_vy((*m_poly).get_vy(i)+(*m_poly).get_fy(i)/(*m_poly).get_poly_mass()*m_integrator_timestep, i);
    (*m_poly_new).set_vz((*m_poly).get_vz(i)+(*m_poly).get_fz(i)/(*m_poly).get_poly_mass()*m_integrator_timestep, i);

  }
   
  m_poly=std::make_unique<Polymer>(*m_poly_new);
}

void Integrator::velocity_verlet()
{

  for(unsigned int i=0; i<(*m_poly).get_poly_sphere(); ++i){
    (*m_poly_new).set_x(pbc((*m_poly).get_x(i)+
                 (*m_poly).get_vx(i)*m_integrator_timestep+
                 0.5*((*m_poly).get_fx(i))/(*m_poly).get_poly_mass()*
                 m_integrator_timestep*m_integrator_timestep), i);

    (*m_poly_new).set_y(pbc((*m_poly).get_y(i)+
                 (*m_poly).get_vy(i)*m_integrator_timestep+
                 0.5*((*m_poly).get_fy(i))/(*m_poly).get_poly_mass()*
                 m_integrator_timestep*m_integrator_timestep), i);
   
    (*m_poly_new).set_z(pbc((*m_poly).get_z(i)+
                 (*m_poly).get_vz(i)*m_integrator_timestep+
                 0.5*((*m_poly).get_fz(i))/(*m_poly).get_poly_mass()*
                 m_integrator_timestep*m_integrator_timestep), i);

    (*m_poly_new).set_vx((*m_poly).get_vx(i)+
                 0.5*((*m_poly).get_fx(i)+(*m_poly_new).get_fx(i))/
                 (*m_poly).get_poly_mass()*m_integrator_timestep, i); 
    
    (*m_poly_new).set_vy((*m_poly).get_vy(i)+
                 0.5*((*m_poly).get_fy(i)+(*m_poly_new).get_fy(i))/
                 (*m_poly).get_poly_mass()*m_integrator_timestep, i); 
    
    (*m_poly_new).set_vz((*m_poly).get_vz(i)+
                 0.5*((*m_poly).get_fz(i)+(*m_poly_new).get_fz(i))/
                 (*m_poly).get_poly_mass()*m_integrator_timestep, i); 
    
  }
    m_poly=std::make_unique<Polymer>(*m_poly_new);
}


void Integrator::langevin_euler()
{
  std::random_device rd;
  std::mt19937 mt (rd());
  std::normal_distribution<double> gauss_term(0., 1.);
  double stoch_term = std::sqrt(2.*m_integrator_gamma*m_integrator_temp/m_integrator_timestep);
  
  for(unsigned int i=0; i<(*m_poly).get_poly_sphere(); ++i){
    (*m_poly_new).set_x(pbc((*m_poly).get_x(i)+
                 (*m_poly).get_vx(i)*m_integrator_timestep), i);
    
    (*m_poly_new).set_y(pbc((*m_poly).get_y(i)+
                 (*m_poly).get_vy(i)*m_integrator_timestep), i);
    
    (*m_poly_new).set_z(pbc((*m_poly).get_z(i)+
                 (*m_poly).get_vz(i)*m_integrator_timestep), i);
    
    (*m_poly_new).set_vx((*m_poly).get_vx(i)+
                 ((*m_poly).get_fx(i)-m_integrator_gamma*(*m_poly).get_vx(i)+stoch_term*gauss_term(mt))/
                 (*m_poly).get_poly_mass()*m_integrator_timestep, i); 
    
    (*m_poly_new).set_vy((*m_poly).get_vy(i)+
                 ((*m_poly).get_fy(i)-m_integrator_gamma*(*m_poly).get_vy(i)+stoch_term*gauss_term(mt))/
                 (*m_poly).get_poly_mass()*m_integrator_timestep, i); 
    
    (*m_poly_new).set_vz((*m_poly).get_vz(i)+
                 ((*m_poly).get_fz(i)-m_integrator_gamma*(*m_poly).get_vz(i)+stoch_term*gauss_term(mt))/
                 (*m_poly).get_poly_mass()*m_integrator_timestep, i); 
    
  }
  m_poly_old=std::make_unique<Polymer>(*m_poly);
  m_poly=std::make_unique<Polymer>(*m_poly_new);
}

void Integrator::langevin_verlet()
{
  std::random_device rd;
  std::mt19937 mt (rd());
  std::normal_distribution<double> gauss_term(0., 1.);
  double stoch_term = std::sqrt(2.*m_integrator_gamma*m_integrator_temp/m_integrator_timestep);
  
  for(unsigned int i=0; i<(*m_poly).get_poly_sphere(); ++i){
    (*m_poly_new).set_x(pbc(2.*(*m_poly).get_x(i)-(*m_poly_old).get_x(i)+
                     ((*m_poly).get_fx(i)-m_integrator_gamma*(*m_poly).get_vx(i)+
                     stoch_term*gauss_term(mt))*m_integrator_timestep*m_integrator_timestep/(*m_poly).get_poly_mass()), i);
    
    (*m_poly_new).set_y(pbc(2.*(*m_poly).get_y(i)-(*m_poly_old).get_y(i)+
                     ((*m_poly).get_fy(i)-m_integrator_gamma*(*m_poly).get_vy(i)+
                     stoch_term*gauss_term(mt))*m_integrator_timestep*m_integrator_timestep/(*m_poly).get_poly_mass()), i);
    
    (*m_poly_new).set_z(pbc(2.*(*m_poly).get_z(i)-(*m_poly_old).get_z(i)+
                     ((*m_poly).get_fz(i)-m_integrator_gamma*(*m_poly).get_vz(i)+
                     stoch_term*gauss_term(mt))*m_integrator_timestep*m_integrator_timestep/(*m_poly).get_poly_mass()), i);
    
    (*m_poly_new).set_vx(pbc((*m_poly_new).get_x(i)-(*m_poly_old).get_x(i))/(2.*m_integrator_timestep), i);
    
    (*m_poly_new).set_vy(pbc((*m_poly_new).get_y(i)-(*m_poly_old).get_y(i))/(2.*m_integrator_timestep), i);
    
    (*m_poly_new).set_vz(pbc((*m_poly_new).get_z(i)-(*m_poly_old).get_z(i))/(2.*m_integrator_timestep), i);
 
  }
  m_poly_old=std::make_unique<Polymer>(*m_poly);
  m_poly=std::make_unique<Polymer>(*m_poly_new);
}


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

    (*m_poly).set_vx(pbc((*m_poly_new).get_x(i)-(*m_poly).get_x(i))/m_integrator_timestep, i); 
    
    (*m_poly).set_vy(pbc((*m_poly_new).get_y(i)-(*m_poly).get_y(i))/m_integrator_timestep, i); 
    
    (*m_poly).set_vz(pbc((*m_poly_new).get_z(i)-(*m_poly).get_z(i))/m_integrator_timestep, i); 
  
    }
  m_poly=std::make_unique<Polymer>(*m_poly_new);
}

void Integrator::scale_factor()
{
  double sum2(0.), scale_factor(0.);
  for( unsigned int i = 0; i<(*m_poly).get_poly_sphere(); ++i){
    
    sum2 += std::pow((*m_poly).get_vx(i),2)+
            std::pow((*m_poly).get_vy(i),2)+
            std::pow((*m_poly).get_vz(i),2);
   
   }
   sum2 /= (double)(*m_poly).get_poly_sphere();
   scale_factor = std::sqrt(3 * m_integrator_temp / (sum2*(*m_poly).get_poly_mass()));   // fs = velocity scale factor 
   for (int i=0; i<(*m_poly).get_poly_sphere(); ++i){
     (*m_poly).set_vx(pbc((*m_poly).get_vx(i) * scale_factor), i);
     (*m_poly).set_vy(pbc((*m_poly).get_vy(i) * scale_factor), i);
     (*m_poly).set_vz(pbc((*m_poly).get_vz(i) * scale_factor), i);
   }
}

