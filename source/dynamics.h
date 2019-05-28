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
#include "extruder.h"

class Dynamics :  public Integrator, public Potential
{
  public: 
    Dynamics(Polymer &poly, std::vector<Extruder> extr, Parameters parm) : 
                                   Integrator(poly, extr, parm), 
                                   Potential(poly, extr,  parm), 
                                   m_dynamics_print(parm.get_print()),
                                   m_dynamics_nstep(parm.get_nstep())  
    { 
      m_parm=parm;
        m_poly = std::make_unique<Polymer>(poly); 
      m_extr.clear();
      	for(auto & i : extr)
	        m_extr.push_back(std::make_unique<Extruder>(i));
    };

    ~Dynamics()
    {
    };

    // Function to get polymer and extruder
    const Polymer & get_poly() const {  return *m_poly; }
    const std::vector<Extruder >  get_extr() const {  
		  std::vector<Extruder> tmp;
		  for (auto & i : m_extr) tmp.push_back(*i);
	    return tmp; 
    }
    Parameters get_parm(){return m_parm;}

    void run();
    void first_fill();
    void fill_extruder();
    bool extr_overlap(Extruder & );

  private:
    Parameters m_parm;
    std::unique_ptr<Polymer> m_poly;
    std::vector<std::unique_ptr<Extruder>> m_extr;
    double k_on = 0.01;
    double k_off = 0.9995;
 // double k_off =1;
    int n_extr = 2;

    int m_dynamics_print = 100;
    int m_dynamics_nstep = 100000;
};
#endif /* DYNAMICS_H_ */
