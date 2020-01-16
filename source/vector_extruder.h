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

#include "field.h"
#include "extruder.h"
#include "parameters.h"
#include "polymer.h"
#include <memory>

struct IdxTime {
    Extruder extr;
    unsigned long int t;
};

class VectorExtruder {

public:

    VectorExtruder() {};

    VectorExtruder(Parameters parm, Extruder &extr, Polymer &poly)
        : m_extr(extr), m_kon(parm.get_kon()), m_koff(parm.get_koff()),
          m_n_max_extr(parm.get_max_extr()), rho0_tot(parm.get_rho0_tot()),
          integrator_timestep(parm.get_timestep()), dist(0., 1.) {
        this->first_fill(poly);
    };

    VectorExtruder(Parameters parm, Extruder &extr, Polymer &poly, FieldAction cohes_poly_int)
        : m_extr(extr), m_kon(parm.get_kon()), m_koff(parm.get_koff()),
          m_n_max_extr(parm.get_max_extr()), rho0_tot(parm.get_rho0_tot()),
          integrator_timestep(parm.get_timestep()), dist(0., 1.) {
        this->first_fill_field(poly, cohes_poly_int);
    };

    // copy constructor
    VectorExtruder(const VectorExtruder &vector_extr)
        : m_extr(vector_extr.m_extr) {
        m_kon = vector_extr.m_kon;
        m_koff = vector_extr.m_koff;
        m_n_max_extr = vector_extr.m_n_max_extr;
        rho0_tot = vector_extr.rho0_tot;
        integrator_timestep = vector_extr.integrator_timestep;
    };

    ~VectorExtruder() {};

    void first_fill(Polymer &);
    void first_fill_field(Polymer &, FieldAction);
    void update(Polymer &);
    void update_diff_density(Polymer &, unsigned long int);
    void update_field(Polymer &, FieldAction);
    double density(Position, unsigned long int, Polymer &);

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

    int vextr_size() {
        return m_vector_extr.size();
    }
private:

    Extruder m_extr;
    Position r, ri;
    std::vector<std::unique_ptr<Extruder>> m_vector_extr;
    std::vector<IdxTime> m_unloaded_extr;

    // maximum number of extruder
    float m_n_max_extr = 3;
    double m_kon = 1E-15;
    double m_koff = 1E-15;

    double integrator_timestep = 1E6;

    // density of extruder in nucleus
    double rho0_tot = 8.9E-2;

    double Nb_eq = 1E2;

    // diffusion of unloaded cohesin from stokes law
    // if cohesin diameter is 0.89a
    const double D = 2.57E-9;
    // move unloaded extruders from local to background
    // every dt step
    int dt = 10000;

    std::mt19937 mt{std::random_device{}()};
    std::uniform_real_distribution<double> dist;
};

#endif /*EXTRUDER_VECTOR_H_*/
