#include "dynamics.h"

void Dynamics::run(bool rouse, bool soft_core, bool lennard_jones, bool compute_energy, std::string output)
{
    for(unsigned long int i=0; i<m_dynamics_nstep; ++i){  
	if(compute_energy) m_poly_old = *m_poly;
        (*m_poly).reset_force();
        if(compute_energy) (*m_poly).reset_energy();

        Potential::set_new_polymer(*m_poly);
        
	if(rouse)
	  this->harmonic_spring_f(compute_energy);
	else if(soft_core){
	  this->harmonic_spring_f(compute_energy);
      	  this->soft_core_f(i, compute_energy);
	}
	else if(lennard_jones){
          this->harmonic_spring_f(compute_energy);
          this->lennard_jones_f(i, compute_energy);
	}
      
        m_poly = std::make_unique<Polymer>(Potential::get_poly());

        Integrator::set_new_polymer(*m_poly);

        this->langevin_overdamped();

        m_poly = std::make_unique<Polymer>(Integrator::get_poly());
        
         if(i%m_dynamics_print == 0) {
            print_xyz(*m_poly, output);
          if(compute_energy) std::cout << i*m_parm.get_timestep()/1E12 << " " << delta_h() << std::endl;
        }
    }
}

void Dynamics::run_extrusion(bool rouse, bool soft_core, bool lennard_jones, bool compute_energy, std::string output) {
  for (unsigned long int i = 0; i < m_dynamics_nstep; ++i) {
    if(compute_energy) m_poly_old = *m_poly;
    if (i % m_dynamics_print == 0)
      m_vector_extr.update(*m_poly);

    (*m_poly).reset_force();
    if(compute_energy) (*m_poly).reset_energy();

    Potential::set_new_polymer(*m_poly);
    Potential::set_new_extruder(m_vector_extr);

    this->extruder_spring_f(compute_energy);
    if(rouse)
      this->harmonic_spring_f(compute_energy);
    else if(soft_core){
      this->harmonic_spring_f(compute_energy);
      this->soft_core_f(i, compute_energy);
    }
    else if(lennard_jones){
      this->harmonic_spring_f(compute_energy);
      this->lennard_jones_f(i, compute_energy);
    }

    m_poly = std::make_unique<Polymer>(Potential::get_poly());
    m_vector_extr = Potential::get_extr();

    Integrator::set_new_polymer(*m_poly);
    Integrator::set_new_extruder(m_vector_extr);

    this->markov_chain();
    this->langevin_overdamped();

    m_poly = std::make_unique<Polymer>(Integrator::get_poly());
    m_vector_extr = Integrator::get_extr();
    
    if (i % m_dynamics_print == 0) {
      print_xyz(*m_poly, output);
    if(compute_energy) std::cout << i*m_parm.get_timestep()/1E12 << " " << delta_h() << std::endl;
      int num_extr = 0;
      for (auto &i : m_vector_extr) {
        print_r(*m_poly, *i,
                 output + std::to_string(num_extr) + ".le");
        ++num_extr;
      }
    }
  }
}

void Dynamics::run_extrusion_field(bool rouse, bool soft_core, bool lennard_jones, bool compute_energy, std::string output) {
  for (unsigned long int i = 0; i < m_dynamics_nstep; ++i) {
    if(compute_energy) m_poly_old = *m_poly;
    if (i % m_dynamics_print == 0)
      m_vector_extr.update_field(*m_poly, m_interaction);
    (*m_poly).reset_force();
    if(compute_energy) (*m_poly).reset_energy();
    Potential::set_new_polymer(*m_poly);
    Potential::set_new_extruder(m_vector_extr);

    this->extruder_spring_f(compute_energy);
    if(rouse)
      this->harmonic_spring_f(compute_energy);
    else if(soft_core){
      this->harmonic_spring_f(compute_energy);
      this->soft_core_f(i, compute_energy);
    }
    else if(lennard_jones){
      this->harmonic_spring_f(compute_energy);
      this->lennard_jones_f(i, compute_energy);
    }

    m_poly = std::make_unique<Polymer>(Potential::get_poly());
    m_vector_extr = Potential::get_extr();
    
    Integrator::set_new_polymer(*m_poly);
    Integrator::set_new_extruder(m_vector_extr);
    Integrator::set_new_field(m_interaction); 

    this->markov_chain();
    this->langevin_overdamped();
    this->extruders_diffusion();
   

    m_poly = std::make_unique<Polymer>(Integrator::get_poly());
    m_vector_extr = Integrator::get_extr();
    m_interaction = Integrator::get_field();
    
    if (i % m_dynamics_print == 0) {
      print_xyz(*m_poly, output);
    if(compute_energy) std::cout << i*m_parm.get_timestep()/1E12 << " " << delta_h() << std::endl;
      int num_extr = 0;
      for (auto &i : m_vector_extr) {
        print_r(*m_poly, *i,
                 output + std::to_string(num_extr) + ".le");
        ++num_extr;
      }
  	print_field(m_interaction, "output/field.xyz");
    }
  }
}

double Dynamics::delta_h(){

    double first_term=0;
    double second_term=0;
    double delta_energy;

    for (unsigned int j=0; j< m_parm.get_nmonomers(); ++j){
      first_term+=m_parm.get_timestep()/(4.*m_parm.get_gamma())*
      ((std::pow((*m_poly).get_fx(j),2)+
       std::pow((*m_poly).get_fy(j),2)+
       std::pow((*m_poly).get_fz(j),2))-
      (std::pow(m_poly_old.get_fx(j),2)+
       std::pow(m_poly_old.get_fy(j),2)+
       std::pow(m_poly_old.get_fz(j),2)));

      second_term+=(((*m_poly).get_x(j)-m_poly_old.get_x(j))*
       ((*m_poly).get_fx(j)+m_poly_old.get_fx(j))/2.+
       ((*m_poly).get_y(j)-m_poly_old.get_y(j))*
       ((*m_poly).get_fy(j)+m_poly_old.get_fy(j))/2.+
       ((*m_poly).get_z(j)-m_poly_old.get_z(j))*
       ((*m_poly).get_fz(j)+m_poly_old.get_fz(j))/2.);
    }

   delta_energy=(*m_poly).get_energy()-m_poly_old.get_energy();

   return first_term+second_term+delta_energy;
}
