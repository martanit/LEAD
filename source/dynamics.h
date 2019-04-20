/*
 * dynamics.h
 *
 *  Created on: March 11, 2019
 *  	Author: martanit
 */

#ifndef DYNAMICS_H_
#define DYNAMICS_H_

#include "polymer.h"

class Dynamics
{
  public: 
    Dynamics(Polymer & poly) : m_poly(poly) { };
    ~Dynamics() { };
    
    void first_move();
    void velocity_verlet();

  private:
    Polymer m_poly;
};

#endif /* DYNAMICS_H_ */
