/*
 * vector_extruder.h
 *
 *  Created on: May 29, 2019
 *  	Author: martanit
 */
#ifndef EXTRUDER_VECTOR_H_
#define EXTRUDER_VECTOR_H_

#include <random>
#include <vector>
#include <memory>

#include "parameters.h"
#include "polymer.h"
#include "extruder.h"

class VectorExtruder {

public:

    VectorExtruder() {};

    VectorExtruder(Parameters parm, Extruder &extr, Polymer &poly)
        : m_extr(extr),
	  m_kon(parm.get_kon()),
	  m_koff(parm.get_koff()),
          integrator_timestep(parm.get_timestep()),
          dist(0., 1.) {
        this->first_fill(poly);
    };

    // copy constructor
    VectorExtruder(const VectorExtruder &vector_extr)
        : m_extr(vector_extr.m_extr) {
        m_kon = vector_extr.m_kon;
        m_koff = vector_extr.m_koff;
        integrator_timestep = vector_extr.integrator_timestep;
    };

    ~VectorExtruder() {};

    // function to fill and update
    // the vector of extrudes
    void first_fill(Polymer &);
    void update(Polymer &);

    // check if two extruders overlap
    friend bool extr_overlap(Extruder &extr);

    bool overlap_lr(Extruder &);
    bool overlap_rl(Extruder &);
    bool overlap_ll(Extruder &);
    bool overlap_rr(Extruder &);

    // define some useful operator
    auto begin() {
        return m_vector_extr.begin();
    }
    auto begin() const {
        return m_vector_extr.begin();
    }
    auto end() {
        return m_vector_extr.end();
    }
    auto end() const {
        return m_vector_extr.end();
    }

    void operator=(const VectorExtruder &lhs) {
        m_vector_extr.clear();
        for (auto &i : (lhs.m_vector_extr))
            m_vector_extr.push_back(std::make_unique<Extruder>(*i));
    }
    
    int vextr_size() {
        return m_vector_extr.size();
    }

private:

    Extruder m_extr;
    std::vector<std::unique_ptr<Extruder>> m_vector_extr;

    double m_kon = 1E-15;
    double m_koff = 1E-15;

    double integrator_timestep = 1E6;

    // random stuff
    std::mt19937 mt{std::random_device{}()};
    std::uniform_real_distribution<double> dist;
};

#endif /*EXTRUDER_VECTOR_H_*/
