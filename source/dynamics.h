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
#include "vector_extruder.h"
#include "integrator.h"
#include "potential.h"

class Dynamics : public Integrator, public Potential {

public:

    Dynamics(Polymer &poly, Parameters parm)
        : Integrator(poly, parm),
	  Potential(poly, parm),
	  m_dynamics_print(parm.get_print()),
          m_dynamics_nstep(parm.get_nstep())
    {
        m_parm = parm;
    };

    Dynamics(Polymer &poly, VectorExtruder &vector_extr, Parameters parm)
        : Integrator(poly, vector_extr, parm),
	  Potential(poly, vector_extr, parm),
	  m_dynamics_print(parm.get_print()),
          m_dynamics_nstep(parm.get_nstep())
    {
        m_parm = parm;
    };

    ~Dynamics() {};

    // run functions
    void run(bool, bool, bool, bool, std::string);
    void run_extrusion(bool, bool, bool, bool, std::string);
    
    // equivalent energy functions
    double delta_h();

    bool print_sys(std::string);

private:

    Parameters m_parm;
    Polymer m_poly_old;
    
    // equivalent energy parameters
    double deltaH, A,B,C;

    // dynamics parameters
    int m_dynamics_print = 100;
    float m_dynamics_nstep = 100000;
};

#endif /* DYNAMICS_H_ */
