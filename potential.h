/*
 * potential.h
 *
 *  Created on: March 11, 2019
 *  	Author: martanit
 */

#ifndef POTENTIAL_H_
#define POTENTIAL_H_

#include "polymer.h"

class Potential
{
  public:
	  Potential(Polymer & poly) : m_poly(poly) { };
	  ~Potential() { };

  private:
	  Polymer m_poly;
};

class LennardJones : public Potential
{
  public:
		LennardJones();
		~LennardJones();
};

class HarmonicSpring : public Potential
{
	public:
		HarmonicSpring();
		~HarmonicSpring();

};

#endif /* POTENTIAL_H_ */
