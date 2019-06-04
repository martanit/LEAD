/*
 * extruder.h
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
        Extruder(Parameters parm, bool rate) : m_parm(parm),
                                    m_rate_vfwl(parm.get_rate_vfwl()),
                                    m_rate_vfwr(parm.get_rate_vfwr()),
                                    m_rate_vbwl(parm.get_rate_vbwl()),
                                    m_rate_vbwr(parm.get_rate_vbwr()),
                                    m_coupling_prob(parm.get_coupling_prob()),
                                    m_ctcf(parm.get_ctcf()),
                                    m_perm_ctcf(parm.get_permeability()),
                                    try_extruder_pos(1, parm.get_psphere()-2),
                                    coupling_try(0,1)
        {
        };
        Extruder(Parameters parm) : m_parm(parm),
                                    m_rate_fwl(parm.get_rate_fwl()),
                                    m_rate_fwr(parm.get_rate_fwr()),
                                    m_rate_bwl(parm.get_rate_bwl()),
                                    m_rate_bwr(parm.get_rate_bwr()),
                                    m_coupling_prob(parm.get_coupling_prob()),
                                    m_ctcf(parm.get_ctcf()),
                                    m_perm_ctcf(parm.get_permeability()),
                                    try_extruder_pos(1, parm.get_psphere()-2),
                                    coupling_try(0,1)
        {
        };

        ~Extruder(){};
        bool operator==(const Extruder& lhs) {
                return lhs.m_extruder_r == this->m_extruder_r and
                       lhs.m_extruder_l == this->m_extruder_l;
            }

        void place_extruder(Polymer poly);
        
        inline const Extruder& get_extr() const{ return *this;};
        inline const int& get_r() const {  return m_extruder_r; };
        inline const int& get_l() const {  return m_extruder_l; };
        inline const std::vector<int>& get_ctcf() const {return m_ctcf;};        
        
        inline const double& get_rate_vfwl(int pos) const {return m_rate_vfwl[pos]; };
        inline const double& get_rate_vfwr(int pos) const {return m_rate_vfwr[pos]; };
        inline const double& get_rate_vbwl(int pos) const {return m_rate_vbwl[pos]; };
        inline const double& get_rate_vbwr(int pos) const {return m_rate_vbwr[pos]; };

        inline const double& get_rate_fwl() const {return m_rate_fwl; };
        inline const double& get_rate_fwr() const {return m_rate_fwr; };
        inline const double& get_rate_bwl() const {return m_rate_bwl; };
        inline const double& get_rate_bwr() const {return m_rate_bwr; };

        
        void set_r(double r){ m_extruder_r = r; } 
        void set_l(double l){ m_extruder_l = l; } 
        
        bool extr_overlap(Extruder & extr);
        
        friend bool print_r(Polymer&, Extruder&, std::string); 
        
    protected:
        int m_extruder_r, m_extruder_l;
        
        std::vector<double> m_rate_vfwl;
        std::vector<double> m_rate_vfwr;
        std::vector<double> m_rate_vbwl;
        std::vector<double> m_rate_vbwr;
     
        double m_rate_fwl = 0.0001;
        double m_rate_fwr = 0.0001;
        double m_rate_bwl = 0.0001;
        double m_rate_bwr = 0.0001;
        double m_perm_ctcf = 0.9;

        std::vector<int> m_ctcf;
        std::vector<double> m_coupling_prob;

    private:
        std::mt19937 mt { std::random_device{}() };    
        std::uniform_int_distribution<> try_extruder_pos;
        std::uniform_real_distribution<double> coupling_try;
        
        Parameters m_parm;
        
        bool set = 0;
        int tmp_extr_pos;
        double tmp_coupling_try;

};  

#endif /*EXTRUDER_H_*/
