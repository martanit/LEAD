#include "dynamics.h"

//integrator

void Dynamics::euler()
{
  for(unsigned int i=0; i<(*m_poly).get_poly_sphere(); ++i){
  
    (*m_poly_new).set_x((*m_poly).get_x(i)+(*m_poly).get_vx(i)*m_timestep, i);
    (*m_poly_new).set_y((*m_poly).get_y(i)+(*m_poly).get_vy(i)*m_timestep, i);
    (*m_poly_new).set_z((*m_poly).get_z(i)+(*m_poly).get_vz(i)*m_timestep, i);
    
    (*m_poly_new).set_vx((*m_poly).get_vx(i)+(*m_poly).get_force_x(i)/(*m_poly).get_poly_mass()*m_timestep, i);
    (*m_poly_new).set_vy((*m_poly).get_vy(i)+(*m_poly).get_force_y(i)/(*m_poly).get_poly_mass()*m_timestep, i);
    (*m_poly_new).set_vz((*m_poly).get_vz(i)+(*m_poly).get_force_z(i)/(*m_poly).get_poly_mass()*m_timestep, i);
  }
   
  m_poly=std::make_unique<Polymer>(*m_poly_new);
}

void Dynamics::velocity_verlet()
{

  for(unsigned int i=0; i<(*m_poly).get_poly_sphere(); ++i){
    (*m_poly_new).set_x(m_conf.pbc((*m_poly).get_x(i)+
                 (*m_poly).get_vx(i)*m_timestep+
                 0.5*((*m_poly).get_force_x(i))/(*m_poly).get_poly_mass()*
                 m_timestep*m_timestep), i);

    (*m_poly_new).set_y(m_conf.pbc((*m_poly).get_y(i)+
                 (*m_poly).get_vy(i)*m_timestep+
                 0.5*((*m_poly).get_force_y(i))/(*m_poly).get_poly_mass()*
                 m_timestep*m_timestep), i);
   
    (*m_poly_new).set_z(m_conf.pbc((*m_poly).get_z(i)+
                 (*m_poly).get_vz(i)*m_timestep+
                 0.5*((*m_poly).get_force_z(i))/(*m_poly).get_poly_mass()*
                 m_timestep*m_timestep), i);

    (*m_poly_new).set_vx((*m_poly).get_vx(i)+
                 0.5*((*m_poly).get_force_x(i)+(*m_poly_new).get_force_x(i))/
                 (*m_poly).get_poly_mass()*m_timestep, i); 
    
    (*m_poly_new).set_vy((*m_poly).get_vy(i)+
                 0.5*((*m_poly).get_force_y(i)+(*m_poly_new).get_force_y(i))/
                 (*m_poly).get_poly_mass()*m_timestep, i); 
    
    (*m_poly_new).set_vz((*m_poly).get_vz(i)+
                 0.5*((*m_poly).get_force_z(i)+(*m_poly_new).get_force_z(i))/
                 (*m_poly).get_poly_mass()*m_timestep, i); 
    
  }
    m_poly=std::make_unique<Polymer>(*m_poly_new);
}


void Dynamics::langevin_euler()
{
  std::random_device rd;
  std::mt19937 mt (rd());
  std::normal_distribution<double> gauss_term(0., 1.);
  double stoch_term = std::sqrt(2.*m_gamma*m_temp/m_timestep);
  
  for(unsigned int i=0; i<(*m_poly).get_poly_sphere(); ++i){
    (*m_poly_new).set_x(m_conf.pbc((*m_poly).get_x(i)+
                 (*m_poly).get_vx(i)*m_timestep), i);
    
    (*m_poly_new).set_y(m_conf.pbc((*m_poly).get_y(i)+
                 (*m_poly).get_vy(i)*m_timestep), i);
    
    (*m_poly_new).set_z(m_conf.pbc((*m_poly).get_z(i)+
                 (*m_poly).get_vz(i)*m_timestep), i);
    
    (*m_poly_new).set_vx((*m_poly).get_vx(i)+
                 ((*m_poly).get_force_x(i)-m_gamma*(*m_poly).get_vx(i)+stoch_term*gauss_term(mt))/
                 (*m_poly).get_poly_mass()*m_timestep, i); 
    
    (*m_poly_new).set_vy((*m_poly).get_vy(i)+
                 ((*m_poly).get_force_y(i)-m_gamma*(*m_poly).get_vy(i)+stoch_term*gauss_term(mt))/
                 (*m_poly).get_poly_mass()*m_timestep, i); 
    
    (*m_poly_new).set_vz((*m_poly).get_vz(i)+
                 ((*m_poly).get_force_z(i)-m_gamma*(*m_poly).get_vz(i)+stoch_term*gauss_term(mt))/
                 (*m_poly).get_poly_mass()*m_timestep, i); 
    
  }
  m_poly_old=std::make_unique<Polymer>(*m_poly);
  m_poly=std::make_unique<Polymer>(*m_poly_new);
}

void Dynamics::langevin_verlet()
{
  std::random_device rd;
  std::mt19937 mt (rd());
  std::normal_distribution<double> gauss_term(0., 1.);
  double stoch_term = std::sqrt(2.*m_gamma*m_temp/m_timestep);
  
  for(unsigned int i=0; i<(*m_poly).get_poly_sphere(); ++i){
    (*m_poly_new).set_x(m_conf.pbc(2.*(*m_poly).get_x(i)-(*m_poly_old).get_x(i)+
                     ((*m_poly).get_force_x(i)-m_gamma*(*m_poly).get_vx(i)+
                     stoch_term*gauss_term(mt))*m_timestep*m_timestep/(*m_poly).get_poly_mass()), i);
    
    (*m_poly_new).set_y(m_conf.pbc(2.*(*m_poly).get_y(i)-(*m_poly_old).get_y(i)+
                     ((*m_poly).get_force_y(i)-m_gamma*(*m_poly).get_vy(i)+
                     stoch_term*gauss_term(mt))*m_timestep*m_timestep/(*m_poly).get_poly_mass()), i);
    
    (*m_poly_new).set_z(m_conf.pbc(2.*(*m_poly).get_z(i)-(*m_poly_old).get_z(i)+
                     ((*m_poly).get_force_z(i)-m_gamma*(*m_poly).get_vz(i)+
                     stoch_term*gauss_term(mt))*m_timestep*m_timestep/(*m_poly).get_poly_mass()), i);
    
    (*m_poly_new).set_vx(m_conf.pbc((*m_poly_new).get_x(i)-(*m_poly_old).get_x(i))/(2.*m_timestep), i);
    
    (*m_poly_new).set_vy(m_conf.pbc((*m_poly_new).get_y(i)-(*m_poly_old).get_y(i))/(2.*m_timestep), i);
    
    (*m_poly_new).set_vz(m_conf.pbc((*m_poly_new).get_z(i)-(*m_poly_old).get_z(i))/(2.*m_timestep), i);
 
  }
  m_poly_old=std::make_unique<Polymer>(*m_poly);
  m_poly=std::make_unique<Polymer>(*m_poly_new);
}


void Dynamics::langevin_overdamped()
{
  std::random_device rd;
  std::mt19937 mt(rd());
  std::normal_distribution<double> gauss_term(0., 1.);
  double stoch_term = std::sqrt(2.*m_gamma*m_temp/m_timestep);

  
  for(unsigned int i=0; i<(*m_poly).get_poly_sphere(); ++i){
    (*m_poly_new).set_x(m_conf.pbc((*m_poly).get_x(i)+
                     ((*m_poly).get_force_x(i)+
                     stoch_term*gauss_term(mt))*m_timestep/m_gamma), i);

    (*m_poly_new).set_y(m_conf.pbc((*m_poly).get_y(i)+
                     ((*m_poly).get_force_y(i)+
                     stoch_term*gauss_term(mt))*m_timestep/m_gamma), i);
   
    (*m_poly_new).set_z(m_conf.pbc((*m_poly).get_z(i)+
                     ((*m_poly).get_force_z(i)+
                     stoch_term*gauss_term(mt))*m_timestep/m_gamma), i);

    (*m_poly).set_vx(m_conf.pbc((*m_poly_new).get_x(i)-(*m_poly).get_x(i))/m_timestep, i); 
    
    (*m_poly).set_vy(m_conf.pbc((*m_poly_new).get_y(i)-(*m_poly).get_y(i))/m_timestep, i); 
    
    (*m_poly).set_vz(m_conf.pbc((*m_poly_new).get_z(i)-(*m_poly).get_z(i))/m_timestep, i); 
  
    }
  m_poly=std::make_unique<Polymer>(*m_poly_new);
}

void Dynamics::scale_factor()
{
  double sum2(0.), scale_factor(0.);
  for( unsigned int i = 0; i<(*m_poly).get_poly_sphere(); ++i){
    
    sum2 += std::pow((*m_poly).get_vx(i),2)+
            std::pow((*m_poly).get_vy(i),2)+
            std::pow((*m_poly).get_vz(i),2);
   
   }
   sum2 /= (double)(*m_poly).get_poly_sphere();
   scale_factor = std::sqrt(3 * m_temp / (sum2*(*m_poly).get_poly_mass()));   // fs = velocity scale factor 
   for (int i=0; i<(*m_poly).get_poly_sphere(); ++i){
     (*m_poly).set_vx(m_conf.pbc((*m_poly).get_vx(i) * scale_factor), i);
     (*m_poly).set_vy(m_conf.pbc((*m_poly).get_vy(i) * scale_factor), i);
     (*m_poly).set_vz(m_conf.pbc((*m_poly).get_vz(i) * scale_factor), i);
   }
}

void Dynamics::run()
{
  for(unsigned int i=0; i<m_parm.get_nstep(); ++i){ 
  //  m_pot.lennard_jones_f();
    m_pot.harmonic_spring_f();
    m_poly = std::make_unique<Polymer>(m_pot.get_poly());
    //if(i==0) this->langevin_euler();
    // else this->langevin_verlet();
    this-> euler();
    m_poly = std::make_unique<Polymer>(this->get_poly());  
    m_pot = Potential(*m_poly, m_parm);
    
    if(i%m_parm.get_print() == 0) print_xyz(*m_poly, "prova.xyz");
    
  }
}
