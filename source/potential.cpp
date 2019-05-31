#include "potential.h"

void Potential::lennard_jones_f()
{
  double x, y, z;
  for (unsigned int i=0; i<m_poly.get_poly_sphere(); ++i )
    for (unsigned int j=i+2; j < m_poly.get_poly_sphere(); ++j){
      x = pbc( m_poly.get_x(i) - m_poly.get_x(j)); 
      y = pbc( m_poly.get_y(i) - m_poly.get_y(j));
      z = pbc( m_poly.get_z(i) - m_poly.get_z(j));
      
      dr = x*x + y*y + z*z;
      
      if(std::sqrt(dr) < m_pot_rcut) {
        lj_x = (x * m_pot_epsilon * (48.0*m_pot_sigma_12/std::pow(dr,7) - 24.0*m_pot_sigma_6/std::pow(dr,4)));
        lj_y = (y * m_pot_epsilon * (48.0*m_pot_sigma_12/std::pow(dr,7) - 24.0*m_pot_sigma_6/std::pow(dr,4)));
        lj_z = (z * m_pot_epsilon * (48.0*m_pot_sigma_12/std::pow(dr,7) - 24.0*m_pot_sigma_6/std::pow(dr,4)));
        
        m_poly.add_force(i, lj_x, lj_y, lj_z);
        m_poly.add_force(j, -lj_x, -lj_y, -lj_z);
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

