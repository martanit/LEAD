#include "potential.h"

void Potential::lennard_jones_f()
{
   forces_local_x.resize(n, 0);
   forces_local_y.resize(n, 0);
   forces_local_z.resize(n, 0);  

  #pragma omp parallel for
    for (unsigned k=0; k<r;k++ ){
      unsigned i = k/n;
      unsigned j = (k%n)+1;
      if(j<=i) {
        i=n-i-2;
        j=n-j-1;
      }
      x = pbc( m_poly.get_x(i) - m_poly.get_x(j)); 
      y = pbc( m_poly.get_y(i) - m_poly.get_y(j));
      z = pbc( m_poly.get_z(i) - m_poly.get_z(j));
      
      dr = x*x + y*y + z*z;
    //  if(dr < (m_pot_rcut*m_pot_rcut)) {
      
        lj_x = (x * m_pot_epsilon * (48.0*m_pot_sigma_12/std::pow(dr,7) - 24.0*m_pot_sigma_6/std::pow(dr,4)));
        lj_y = (y * m_pot_epsilon * (48.0*m_pot_sigma_12/std::pow(dr,7) - 24.0*m_pot_sigma_6/std::pow(dr,4)));
        lj_z = (z * m_pot_epsilon * (48.0*m_pot_sigma_12/std::pow(dr,7) - 24.0*m_pot_sigma_6/std::pow(dr,4)));

        forces_local_x[i] += lj_x;
        forces_local_y[i] += lj_y;
        forces_local_z[i] += lj_z;
        
        forces_local_x[j] -= lj_x;
        forces_local_y[j] -= lj_y;
        forces_local_z[j] -= lj_z;
     // }
    }
    for(unsigned i=0; i<n; i++)
      m_poly.add_force(i, forces_local_x[i], forces_local_y[i], forces_local_z[i]);
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
  //  for(unsigned i=0; i<10; i++)
   //   std::cout<< m_poly.get_fx(i) << " " << m_poly.get_fy(i) << " " << m_poly.get_fz(i) << std::endl;
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

