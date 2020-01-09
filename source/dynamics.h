/*
 * dynamics.h
 *
 *  Created on: March 11, 2019
 *  	Author: martanit
 */

#ifndef DYNAMICS_H_
#define DYNAMICS_H_

#include "extruder.h"
#include "integrator.h"
#include "parameters.h"
#include "polymer.h"
#include "potential.h"
#include "vector_extruder.h"

class Dynamics : public Integrator, public Potential {
public:
    Dynamics(Polymer &poly, Parameters parm)
        : Integrator(poly, parm), Potential(poly, parm),
          m_dynamics_print(parm.get_print()),
          m_dynamics_nstep(parm.get_nstep())

    {
        m_parm = parm;
        m_poly = std::make_unique<Polymer>(poly);
    };
    Dynamics(Polymer &poly, VectorExtruder &vector_extr, Parameters parm)
        : Integrator(poly, vector_extr, parm), Potential(poly, vector_extr, parm),
          m_vector_extr(vector_extr), m_dynamics_print(parm.get_print()),
          m_dynamics_nstep(parm.get_nstep())

    {
        m_parm = parm;
        m_poly = std::make_unique<Polymer>(poly);
    };

    ~Dynamics() {};

    // Function to get polymer and extruder
    const Polymer &get_poly() const {
        return *m_poly;
    }
    const VectorExtruder &get_extr() const {
        return m_vector_extr;
    }

    void run(bool, bool, bool, bool, std::string);
    void run_extrusion(bool, bool, bool, bool, std::string);
    double delta_h();


private:
    Parameters m_parm;
    std::unique_ptr<Polymer> m_poly;
    VectorExtruder m_vector_extr;
    Polymer m_poly_old;
    double deltaH, A,B,C;

    int m_dynamics_print = 100;
    float m_dynamics_nstep = 100000;
};
#endif /* DYNAMICS_H_ */
