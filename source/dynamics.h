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
    Dynamics(Polymer poly, Dynamics_Parameters dyn_parm, Potential pot) : 
        m_conf(dyn_parm),  
        m_pot(pot),
        m_dyn_timestep(dyn_parm.get_timestep()),
        m_dyn_print(dyn_parm.get_print()),
        m_dyn_nstep(dyn_parm.get_nstep()),
        m_dyn_gamma(dyn_parm.get_gamma()),
        m_dyn_temp(dyn_parm.get_temp())
    {
      m_poly = std::make_unique<Polymer>(poly);
      m_poly_new = std::make_unique<Polymer>(*m_poly);
      m_poly_old = std::make_unique<Polymer>(*m_poly);
    };

    ~Dynamics()
    {
    //  if(m_poly != nullptr) delete m_poly;
    //  if(m_poly_new != nullptr) delete m_poly_new;
    //  if(m_poly_old != nullptr) delete m_poly_old;
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
    
    double m_dyn_timestep = 1.;
    int m_dyn_print = 100;
    int m_dyn_nstep = 100000;
    double m_dyn_gamma=1.;
    double m_dyn_temp = 1.;
};
#endif /* DYNAMICS_H_ */
