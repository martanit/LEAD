/*
 * dynamics.h
 *
 *  Created on: March 11, 2019
 *  	Author: martanit
 */

#ifndef DYNAMICS_H_
#define DYNAMICS_H_

#include "parameters.h"
#include "polymer.h"
#include "potential.h"
#include "integrator.h"

class Dynamics :  public Integrator, public Potential
{
  public: 
    Dynamics(Polymer &poly, Parameters parm) : 
                                   Integrator(poly, parm), 
                                   Potential(poly, parm), 
                                   m_dynamics_print(parm.get_print()),
                                   m_dynamics_nstep(parm.get_nstep())  
    {
        m_poly = std::make_unique<Polymer>(poly); 
    };

    ~Dynamics()
    {
    };

    // Function to get polymer
    const Polymer & get_poly() const {  return *m_poly; }
    void run();

  private:

    std::unique_ptr<Polymer> m_poly;
    
    int m_dynamics_print = 100;
    int m_dynamics_nstep = 100000;
};
#endif /* DYNAMICS_H_ */
