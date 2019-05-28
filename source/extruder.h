/*
 * dynamics.h
 *
 *  Created on: May 20, 2019
 *  	Author: martanit
 */
#ifndef EXTRUDER_H_
#define EXTRUDER_H_

#include <random>
#include <vector>
#include <fstream>

#include "parameters.h"
#include "polymer.h"
#include<memory>

class Extruder  {

    public:
        Extruder(Parameters parm) : m_parm(parm),
                                    m_rate_l(parm.get_rate_l()),
                                    m_rate_r(parm.get_rate_r()),
                                    m_coupling_prob(parm.get_coupling_prob()),
                                    m_ctcf(parm.get_ctcf())
        {
        };

        Extruder(Parameters parm, Polymer poly) : m_parm(parm),
                                    m_rate_l(parm.get_rate_l()),
                                    m_rate_r(parm.get_rate_r()),
                                    m_coupling_prob(parm.get_coupling_prob()),
                                    m_ctcf(parm.get_ctcf())
        {
            this -> place_extruder(poly);
        };

        ~Extruder(){};
        void place_extruder(Polymer poly);
        
        const int& get_r() const {  return m_extruder_r; };
        const int& get_l() const {  return m_extruder_l; };
        
        const double& get_rate_l(int pos) const {return m_rate_l[pos]; };
        const double& get_rate_r(int pos) const {return m_rate_r[pos]; };

        const std::vector<int>& get_ctcf() const {return m_ctcf;};        
        void set_r(double r){ m_extruder_r = r; } 
        void set_l(double l){ m_extruder_l = l; } 
        
        bool extr_overlap(Extruder & extr);
        bool extr_overlap_l(std::vector<std::unique_ptr<Extruder>> & extr);
        bool extr_overlap_r(std::vector<std::unique_ptr<Extruder>> & extr);
        friend bool print_r(Polymer&, Extruder&, std::string); 
    protected:
        int m_extruder_r, m_extruder_l;
        
        std::vector<double> m_rate_l;
        std::vector<double> m_rate_r;
        
        std::vector<int> m_ctcf;
        std::vector<double> m_coupling_prob;
        Parameters m_parm;

};  

#endif /*EXTRUDER_H_*/
