/*
 * dynamics.h
 *
 *  Created on: March 11, 2019
 *  	Author: martanit
 */

#ifndef DYNAMICS_H_
#define DYNAMICS_H_

#include "polymer.h"
#include "potential.h"
#include <cmath>
#include <memory>

class Dynamics
{
  public: 
    Dynamics(Polymer poly, Parameters parm) :  
                                                m_conf(parm),       
                                                m_pot(poly, parm),
                                                m_timestep(parm.get_timestep()),
                                                m_temp(parm.get_temp()),
                                                m_gamma(parm.get_gamma())
    {
     m_parm=parm;
      m_poly = std::make_unique<Polymer>(poly);
      m_poly_new = std::make_unique<Polymer>(*m_poly);
      m_poly_old = std::make_unique<Polymer>(*m_poly);
    };

    ~Dynamics()
    {
  //    if(m_poly != nullptr) delete m_poly;
 //     if(m_poly_new != nullptr) delete m_poly_new;
//      if(m_poly_old != nullptr) delete m_poly_old;
    };

    // Function to get polymer
    const Polymer & get_poly() const {  return *m_poly; }
    
    // Integrators
    void euler();
    void velocity_verlet();
    void langevin_euler();
    void langevin_verlet();
    void langevin_overdamped();
    
    // Thermostat
    void scale_factor();

    // Run simulation
    void run();

  private:
    std::unique_ptr<Polymer> m_poly;
    std::unique_ptr<Polymer> m_poly_new;
    std::unique_ptr<Polymer> m_poly_old;

    Utils m_conf;

    Potential m_pot;
    Parameters m_parm;

    float m_timestep = 1.;
    float m_temp = 1.;
    double m_rcut=2.5;
    double m_gamma=1.;
};
#endif /* DYNAMICS_H_ */
