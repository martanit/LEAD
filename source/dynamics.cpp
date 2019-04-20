#include "dynamics.h"

void Dynamics::velocity_verlet()
{
  for(unsigned int i=0; i<m_poly.m_parm.get_psphere(); ++i){
    m_poly(i,0) = m_poly(i,0)+m_poly(i,3)*m_poly.m_parm.get_timestep()+0.5*a(i)*m_poly.m_parm.get_timestep()*m_poly.m_parm.get_timestep();
    m_poly(i,1) = m_poly(i,1)+m_poly(i,4)*m_poly.m_parm.get_timestep()+0.5*a(i)*m_poly.m_parm.get_timestep()*m_poly.m_parm.get_timestep();
    m_poly(i,2) = m_poly(i,2)+m_poly(i,5)*m_poly.m_parm.get_timestep()+0.5*a(i)*m_poly.m_parm.get_timestep()*m_poly.m_parm.get_timestep();
  
  }
}
