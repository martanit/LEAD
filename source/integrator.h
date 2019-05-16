/*
 * integrator.h
 *
 *  Created on: March 11, 2019
 *  	Author: martanit
 */

#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

#include "polymer.h"
#include "potential.h"
#include <cmath>
#include <memory>

class Integrator
{
  public: 
    Integrator(Polymer poly, Parameters integrator_parm) : 
        m_conf(integrator_parm),  
        m_integrator_timestep(integrator_parm.get_timestep()),
        m_integrator_gamma(integrator_parm.get_gamma()),
        m_integrator_temp(integrator_parm.get_temp())
    {
      m_poly = std::make_unique<Polymer>(poly);
      m_poly_new = std::make_unique<Polymer>(*m_poly);
      m_poly_old = std::make_unique<Polymer>(*m_poly);
    };

    ~Integrator()
    {
    //  if(m_poly != nullptr) delete m_poly;
    //  if(m_poly_new != nullptr) delete m_poly_new;
    //  if(m_poly_old != nullptr) delete m_poly_old;
    };

    // Function to get polymer
    void set_new_polymer(Polymer& poly) { *m_poly = poly;}
    const Polymer & get_poly() const {  return *m_poly; }
    
    // Integrators
    void euler();
    void velocity_verlet();
    void langevin_euler();
    void langevin_verlet();
    void langevin_overdamped();
    
    // Thermostat
    void scale_factor();
  
  private:
    std::unique_ptr<Polymer> m_poly;
    std::unique_ptr<Polymer> m_poly_new;
    std::unique_ptr<Polymer> m_poly_old;

    Utils m_conf;
    
    double m_integrator_timestep = 0.0001;
    double m_integrator_gamma=1.;
    double m_integrator_temp = 1.;
};
#endif /* INTEGRATOR_H_ */
