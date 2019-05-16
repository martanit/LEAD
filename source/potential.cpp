#include "potential.h"

void Potential::lennard_jones_f()
{
  double fx(0.),
         fy(0.),
         fz(0.);
//  m_poly.reset_force();
  double dr;
  std::vector<double> r(3, 0.0);
  for (unsigned int i=0; i<m_poly.get_poly_sphere(); ++i )
    for (unsigned int j=i+2; j < m_poly.get_poly_sphere(); ++j){
      r[0] = m_conf.pbc( m_poly.get_x(i) - m_poly.get_x(j) ); 
      r[1] = m_conf.pbc( m_poly.get_y(i) - m_poly.get_y(j) );
      r[2] = m_conf.pbc( m_poly.get_z(i) - m_poly.get_z(j) );
      
      dr = std::sqrt( std::pow(r[0],2) + std::pow(r[1],2) + std::pow(r[2],2));
      
      if(dr < m_pot_rcut) {
        fx = (r[0] * (48.0/std::pow(dr,14) - 24.0/std::pow(dr,8)));
        fy = (r[1] * (48.0/std::pow(dr,14) - 24.0/std::pow(dr,8)));
        fz = (r[2] * (48.0/std::pow(dr,14) - 24.0/std::pow(dr,8)));
        
        m_poly.add_force(i, fx, fy, fz);
        m_poly.add_force(j, -fx, -fy, -fz);
      }
    }
}

void  Potential::harmonic_spring_f()
{
  double spring_x(0.),
         spring_y(0.),
         spring_z(0.);
//  m_poly.reset_force();
  for (int i=0; i<m_poly.get_poly_sphere()-1; ++i ){
      spring_x = -k * (m_conf.pbc( m_poly.dist(i, i+1) - m_poly.get_poly_dist())) * 
              ( m_conf.pbc( m_poly.get_x(i) - m_poly.get_x(i+1)))/m_poly.dist(i, i+1);
      spring_y = -k * (m_conf.pbc( m_poly.dist(i, i+1) - m_poly.get_poly_dist())) * 
              ( m_conf.pbc( m_poly.get_y(i) - m_poly.get_y(i+1)))/m_poly.dist(i, i+1);
      spring_z = -k * (m_conf.pbc( m_poly.dist(i, i+1) - m_poly.get_poly_dist())) * 
              ( m_conf.pbc( m_poly.get_z(i) - m_poly.get_z(i+1)))/m_poly.dist(i, i+1);
      
      m_poly.add_force(i,spring_x, spring_y, spring_z);
      m_poly.add_force(i+1, -spring_x, -spring_y, -spring_z);
  }
}

