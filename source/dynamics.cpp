#include "dynamics.h"

void Dynamics::run()
{
  Polymer m_poly_old;
  double deltaH, A, B, C;
    for(unsigned int i=0; i<m_dynamics_nstep; ++i){ 
       // if(i%m_dynamics_print == 0)    
       //     m_vector_extr.update(*m_poly);

        m_poly_old = *m_poly;
        (*m_poly).reset_force();
        (*m_poly).reset_energy();
       
        Potential::set_new_polymer(*m_poly);
        //Potential::set_new_extruder(m_vector_extr);
      
      //  this->extruder_spring_f();
        this->lennard_jones_f(i, true);
        this->harmonic_spring_f();
        this->kinetic();
        
        
        m_poly = std::make_unique<Polymer>(Potential::get_poly());
       // m_vector_extr = Potential::get_extr();

        Integrator::set_new_polymer(*m_poly);
      //  Integrator::set_new_extruder(m_vector_extr);
      
      //  this->markov_chain();
        this->langevin_overdamped();

        m_poly = std::make_unique<Polymer>(Integrator::get_poly()); 
     //   m_vector_extr = Integrator::get_extr();
    A=0;
    B=0; 
      for (unsigned int j=0; j< m_parm.get_psphere(); ++j){  
          A+=m_parm.get_timestep()/(4.*m_parm.get_gamma())*
           ((std::pow((*m_poly).get_fx(j),2)+
            std::pow((*m_poly).get_fy(j),2)+
            std::pow((*m_poly).get_fz(j),2))-
           (std::pow(m_poly_old.get_fx(j),2)+
            std::pow(m_poly_old.get_fy(j),2)+
            std::pow(m_poly_old.get_fz(j),2)));
            
           B+=(((*m_poly).get_x(j)-m_poly_old.get_x(j))*
            ((*m_poly).get_fx(j)+m_poly_old.get_fx(j))/2.+
            ((*m_poly).get_y(j)-m_poly_old.get_y(j))*
            ((*m_poly).get_fy(j)+m_poly_old.get_fy(j))/2.+
            ((*m_poly).get_z(j)-m_poly_old.get_z(j))*
            ((*m_poly).get_fz(j)+m_poly_old.get_fz(j))/2.);
      }
      
        C=(*m_poly).get_energy()-m_poly_old.get_energy();
        deltaH=A+B+C;
        std::cout  << i <<" "<< A<<" "<< B <<" "<< C <<" "<< deltaH << std::endl;
        if(i%m_dynamics_print == 0) {
            print_xyz(*m_poly, "output/traj.xyz");
         //   int num_extr=0;
           // for (auto &i : m_vector_extr){
           //   print_r(*m_poly, *i, "output/loop_extrusion_"+std::to_string(num_extr)+".r");
           //   ++num_extr;
           // }
        }
    }
}

