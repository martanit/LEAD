/*
 * integrator.h
 *
 *  Created on: March 11, 2019
 *  	Author: martanit
 */

#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

#include <cmath>
#include <memory>

#include "polymer.h"
#include "potential.h"
#include "extruder.h"

class Integrator
{
  public: 
    Integrator(Polymer poly, Extruder extr, Parameters integrator_parm) : 
        m_integrator_timestep(integrator_parm.get_timestep()),
        m_integrator_gamma(integrator_parm.get_gamma()),
        m_integrator_temp(integrator_parm.get_temp()),
        m_box(integrator_parm.get_box())
    {
      m_poly = std::make_unique<Polymer>(poly);
      m_poly_new = std::make_unique<Polymer>(*m_poly);
      m_poly_old = std::make_unique<Polymer>(*m_poly);
      m_extr = std::make_unique<Extruder>(extr);
    };

    ~Integrator()
    {
    /*    delete m_poly.release();
        delete m_poly_new.release();
        delete m_poly_old.release();
        delete m_extr.release();
   */ };

    // Function to get polymer
<<<<<<< HEAD
    inline void set_new_polymer(Polymer& poly) { *m_poly = poly;}
    inline const Polymer & get_poly() const {  return *m_poly; }
=======
    void set_new_polymer(Polymer& poly) { m_poly = std::make_unique<Polymer>(poly);}
    void set_new_extruder(Extruder& extr) { m_extr = std::make_unique<Extruder>(extr);}
    const Polymer & get_poly() const {  return *m_poly; }
    const Extruder & get_extr() const {  return *m_extr; }
>>>>>>> 608e69bd3f7f70a81d7400f213e161aff156fe32
    
    // polymer integrators
    void euler();
    void velocity_verlet();
    void langevin_euler();
    void langevin_verlet();
    void langevin_overdamped();
    
    // extruder move
    void markov_chain();

    // Thermostat
    void scale_factor();
    
  protected:
    std::unique_ptr<Polymer> m_poly;
    std::unique_ptr<Polymer> m_poly_new;
    std::unique_ptr<Polymer> m_poly_old;
    
    std::unique_ptr<Extruder> m_extr;

    double m_integrator_timestep = 0.0001;
    double m_integrator_gamma=1.;
    double m_integrator_temp = 1.;
    
    double m_box = 50.;
};
#endif /* INTEGRATOR_H_ */
