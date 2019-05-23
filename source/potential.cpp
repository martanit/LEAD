#include "potential.h"

void Potential::lennard_jones_f()
{
  double lj_x(0.),
         lj_y(0.),
         lj_z(0.);
  
  double dr;
  std::vector<double> r(3, 0.0);
  for (unsigned int i=0; i<m_poly.get_poly_sphere(); ++i )
    for (unsigned int j=i+2; j < m_poly.get_poly_sphere(); ++j){
      r[0] = pbc( m_poly.get_x(i) - m_poly.get_x(j)); 
      r[1] = pbc( m_poly.get_y(i) - m_poly.get_y(j));
      r[2] = pbc( m_poly.get_z(i) - m_poly.get_z(j));
      
      dr = std::sqrt( std::pow(r[0],2) + std::pow(r[1],2) + std::pow(r[2],2));
      
      if(dr < m_pot_rcut) {
        lj_x = (r[0] * m_pot_epsilon * (48.0*std::pow(m_pot_sigma,12)/std::pow(dr,14) - 24.0*std::pow(m_pot_sigma,6)/std::pow(dr,8)));
        lj_y = (r[1] * m_pot_epsilon * (48.0*std::pow(m_pot_sigma,12)/std::pow(dr,14) - 24.0*std::pow(m_pot_sigma,6)/std::pow(dr,8)));
        lj_z = (r[2] * m_pot_epsilon * (48.0*std::pow(m_pot_sigma,12)/std::pow(dr,14) - 24.0*std::pow(m_pot_sigma,6)/std::pow(dr,8)));
        
        m_poly.add_force(i, lj_x, lj_y, lj_z);
        m_poly.add_force(j, -lj_x, -lj_y, -lj_z);
      }
    }
}

void  Potential::harmonic_spring_f()
{
  double spring_x(0.),
         spring_y(0.),
         spring_z(0.);
  
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
  double spring_x(0.),
         spring_y(0.),
         spring_z(0.);
      spring_x = -k_extr * pbc( m_poly.dist(m_extr.get_l(), m_extr.get_r()) - extr_lenght) * 
                      pbc( m_poly.get_x(m_extr.get_l()) - m_poly.get_x(m_extr.get_r()))
                      /m_poly.dist(m_extr.get_l(), m_extr.get_r());
      spring_y = -k_extr * pbc( m_poly.dist(m_extr.get_l(), m_extr.get_r()) - extr_lenght) * 
                      pbc( m_poly.get_y(m_extr.get_l()) - m_poly.get_y(m_extr.get_r()))
                      /m_poly.dist(m_extr.get_l(), m_extr.get_r());
      spring_z = -k_extr * pbc( m_poly.dist(m_extr.get_l(), m_extr.get_r()) - extr_lenght) * 
                      pbc( m_poly.get_z(m_extr.get_l()) - m_poly.get_z(m_extr.get_r()))
                      /m_poly.dist(m_extr.get_l(), m_extr.get_r());
      
      m_poly.add_force(m_extr.get_l(),spring_x, spring_y, spring_z);
      m_poly.add_force(m_extr.get_r(), -spring_x, -spring_y, -spring_z);
}

