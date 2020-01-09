/*
 * vector_extruder.h
 *
 *  Created on: May 29, 2019
 *  	Author: martanit
 */
#ifndef EXTRUDER_VECTOR_H_
#define EXTRUDER_VECTOR_H_

#include <fstream>
#include <random>
#include <vector>

#include "extruder.h"
#include "parameters.h"
#include "polymer.h"
#include <memory>

class VectorExtruder {

public:

    VectorExtruder() {};

    VectorExtruder(Parameters parm, Extruder &extr, Polymer &poly)
        : m_extr(extr), m_kon(parm.get_kon()), m_koff(parm.get_koff()),
          m_n_max_extr(parm.get_max_extr()),
          integrator_timestep(parm.get_timestep()), dist(0., 1.) {
        this->first_fill(poly);
    };

    // copy constructor
    VectorExtruder(const VectorExtruder &vector_extr)
        : m_extr(vector_extr.m_extr) {
        m_kon = vector_extr.m_kon;
        m_koff = vector_extr.m_koff;
        m_n_max_extr = vector_extr.m_n_max_extr;
        integrator_timestep = vector_extr.integrator_timestep;
    };

    ~VectorExtruder() {};

    void first_fill(Polymer &);
    void update(Polymer &);

    friend bool extr_overlap(Extruder &extr);

    bool overlap_lr(Extruder &);
    bool overlap_rl(Extruder &);
    bool overlap_ll(Extruder &);
    bool overlap_rr(Extruder &);

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

private:

    Extruder m_extr;
    std::vector<std::unique_ptr<Extruder>> m_vector_extr;

    // maximum number of extruder
    float m_n_max_extr = 3;
    double m_kon = 0.9;
    double m_koff = 0.001;
    double integrator_timestep = 0;

    std::mt19937 mt{std::random_device{}()};
    std::uniform_real_distribution<double> dist;
};

#endif /*EXTRUDER_VECTOR_H_*/
