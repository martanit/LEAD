#include "potential.h"

void Potential::kinetic()
{
    double k_x, k_y, k_z;
    for(unsigned int i=0; i< m_poly.get_poly_sphere(); ++i){

        k_x = m_poly.get_poly_mass()*m_poly.get_vx(i)*m_poly.get_vx(i)/2.;
        k_y = m_poly.get_poly_mass()*m_poly.get_vy(i)*m_poly.get_vy(i)/2.;
        k_z = m_poly.get_poly_mass()*m_poly.get_vz(i)*m_poly.get_vz(i)/2.;
        
        m_poly.add_energy(k_x+k_y+k_z);
    }
}
void Potential::lennard_jones_f(int step, bool attarctive)
{
  
  if(step%100==0 or step==0){
    for (unsigned int i=0; i<m_poly.get_poly_sphere(); ++i ){
        sphere[i].clear();
        for (unsigned int j=i+2; j < m_poly.get_poly_sphere(); ++j){
            
            x = pbc( m_poly.get_x(i) - m_poly.get_x(j)); 
            y = pbc( m_poly.get_y(i) - m_poly.get_y(j));
            z = pbc( m_poly.get_z(i) - m_poly.get_z(j));
            
            dr = x*x + y*y + z*z;

            if(dr < m_pot_rcut*m_pot_rcut){
                f_x = x * m_pot_epsilon * 48.0*m_pot_sigma_12/std::pow(dr,7);
                f_y = y * m_pot_epsilon * 48.0*m_pot_sigma_12/std::pow(dr,7);
                f_z = z * m_pot_epsilon * 48.0*m_pot_sigma_12/std::pow(dr,7);
                
                if(attractive){
                    f_x -=  x * m_pot_epsilon*24*m_pot_sigma_6/std::pow(dr,4);
                    f_y -=  y * m_pot_epsilon*24*m_pot_sigma_6/std::pow(dr,4);
                    f_z -=  z * m_pot_epsilon*24*m_pot_sigma_6/std::pow(dr,4);
                   
                    m_poly.add_energy(-4*m_pot_epsilon*m_pot_sigma_6/std::pow(dr,3));
                }

                m_poly.add_force(i, f_x, f_y, f_z);
                m_poly.add_force(k, -f_x, -f_y, -f_z);
                
                m_poly.add_energy(4*m_pot_epsilon*m_pot_sigma_12/std::pow(dr,6));
                // calculate points that are into m_pot_rcut^2
                sphere[i].push_back(j);
            }
        }
    }
  }
  else{ 
    for(unsigned int i=0; i<m_poly.get_poly_sphere(); ++i){
            for (auto && k : sphere[i]){
                
                x = pbc( m_poly.get_x(i) - m_poly.get_x(k)); 
                y = pbc( m_poly.get_y(i) - m_poly.get_y(k));
                z = pbc( m_poly.get_z(i) - m_poly.get_z(k));
            
                dr = x*x + y*y + z*z;
                if(dr < m_pot_rcut*m_pot_rcut){
                    f_x = x * m_pot_epsilon * 48.0*m_pot_sigma_12/std::pow(dr,7);
                    f_y = y * m_pot_epsilon * 48.0*m_pot_sigma_12/std::pow(dr,7); 
                    f_z = z * m_pot_epsilon * 48.0*m_pot_sigma_12/std::pow(dr,7);
                
                    if(attractive){
                        f_x -=  x * m_pot_epsilon*24*m_pot_sigma_6/std::pow(dr,4);
                        f_y -=  y * m_pot_epsilon*24*m_pot_sigma_6/std::pow(dr,4);
                        f_z -=  z * m_pot_epsilon*24*m_pot_sigma_6/std::pow(dr,4);
                        
                        m_poly.add_energy(-4*m_pot_epsilon*m_pot_sigma_6/std::pow(dr,3));
                    }
                   
                    m_poly.add_force(i, f_x, f_y, f_z);
                    m_poly.add_force(k, -f_x, -f_y, -f_z);
                    
                    m_poly.add_energy(4*m_pot_epsilon*m_pot_sigma_12/std::pow(dr,6));
                }
            }
        }
    }
}

void  Potential::harmonic_spring_f()
{
  for (int i=0; i<m_poly.get_poly_sphere()-1; ++i ){
      spring_x = -k * pbc( m_poly.dist(i, i+1) - m_pot_sigma) * 
                      pbc( m_poly.get_x(i) - m_poly.get_x(i+1))
                      /m_poly.dist(i, i+1);
      spring_y = -k * pbc( m_poly.dist(i, i+1) - m_pot_sigma) * 
                      pbc( m_poly.get_y(i) - m_poly.get_y(i+1))
                      /m_poly.dist(i, i+1);
      spring_z = -k * pbc( m_poly.dist(i, i+1) - m_pot_sigma) * 
                      pbc( m_poly.get_z(i) - m_poly.get_z(i+1))
                      /m_poly.dist(i, i+1);
      
      m_poly.add_force(i,spring_x, spring_y, spring_z);
      m_poly.add_force(i+1, -spring_x, -spring_y, -spring_z);

      m_poly.add_energy(k*pbc( m_poly.dist(i, i+1) - m_pot_sigma)*
                          pbc( m_poly.dist(i, i+1) - m_pot_sigma)/2.); 
  }
}

void Potential::extruder_spring_f()
{  
  //iterate over each extruder
     for(auto &i: m_vector_extr) {
        spring_x = -k_extr * pbc( m_poly.dist((*i).get_l(), (*i).get_r()) - extr_lenght) * 
                      pbc( m_poly.get_x((*i).get_l()) - m_poly.get_x((*i).get_r()))
                      /m_poly.dist((*i).get_l(), (*i).get_r());
        spring_y = -k_extr * pbc( m_poly.dist((*i).get_l(), (*i).get_r()) - extr_lenght) * 
                      pbc( m_poly.get_y((*i).get_l()) - m_poly.get_y((*i).get_r()))
                      /m_poly.dist((*i).get_l(), (*i).get_r());
        spring_z = -k_extr * pbc( m_poly.dist((*i).get_l(), (*i).get_r()) - extr_lenght) * 
                      pbc( m_poly.get_z((*i).get_l()) - m_poly.get_z((*i).get_r()))
                      /m_poly.dist((*i).get_l(), (*i).get_r());
      
        m_poly.add_force((*i).get_l(),spring_x, spring_y, spring_z);
        m_poly.add_force((*i).get_r(), -spring_x, -spring_y, -spring_z);
    }
}






