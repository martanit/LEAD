#include "dynamics.h"

void Dynamics::run()
{
    

    for(unsigned int i=0; i<m_dynamics_nstep; ++i){ 
        if(i%m_dynamics_print == 0){   
          
            this->fill_extruder();
        } 
        
        (*m_poly).reset_force();
       
        Potential::set_new_polymer(*m_poly);
        Potential::set_new_extruder(m_extr);
      
        this->extruder_spring_f();
        this->lennard_jones_f();
        this->harmonic_spring_f();
       
        m_poly = std::make_unique<Polymer>(Potential::get_poly());
        
        m_extr.clear(); 
	      for (auto &i : Potential::get_extr())
		      m_extr.push_back(std::make_unique<Extruder>(i)); 

        Integrator::set_new_polymer(*m_poly);
        Integrator::set_new_extruder(m_extr);
      
        this->markov_chain();
        this->langevin_overdamped();
        
        m_poly = std::make_unique<Polymer>(Integrator::get_poly()); 
        
        m_extr.clear(); 
	      for (const auto &i : Integrator::get_extr())
		      m_extr.push_back(std::make_unique<Extruder>(i)); 

        if(i%m_dynamics_print == 0) {
            print_xyz(*m_poly, "output/traj.xyz");
            int num_extr=0;
            for (auto &i : m_extr){
              print_r(*m_poly, *i, "output/loop_extrusion_"+std::to_string(num_extr)+".r");
              ++num_extr;
              }
          }
    }
}

void Dynamics::fill_extruder()
{
    std::random_device rd;
    std::mt19937 mt (rd());
    std::uniform_real_distribution<double> dist(0., 1.);
    int n_extr=2;
    bool is_overl=false;
    std::vector<Extruder> tmp_extruder;
    for(int i = 0; i<n_extr; ++i){
        Extruder e(this->get_parm(), *m_poly); 
        for(auto & j : m_extr)
            if((*j).extr_overlap(e)) {
                is_overl=true;
                break;
            }
        if( k_on > dist(mt) and !(is_overl)) 
            tmp_extruder.push_back(e);
    }
    m_extr.clear();
    for( auto & i : tmp_extruder)
        m_extr.push_back(std::make_unique<Extruder>(i));
}
