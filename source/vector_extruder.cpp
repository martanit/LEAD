#include "vector_extruder.h"

void VectorExtruder::first_fill(Polymer &poly) {
    bool is_overl = false;
    if(m_extr.place_extruder(poly)){
        if(m_vector_extr.size() != 0) 
            for (auto &j : m_vector_extr)
                if ((*j).extr_overlap(m_extr)) {
                    is_overl = true;
                    break;
                }
        // tmp 10000
        if (poly.get_poly_nmonomers() * m_kon * 0.032 
			* integrator_timestep * 10000  > dist(mt) and !(is_overl))
            m_vector_extr.push_back(std::make_unique<Extruder>(m_extr));
    }
}

void VectorExtruder::update(Polymer &poly) {
    if(m_vector_extr.size() == 0) this->first_fill(poly);
    else {
        bool is_overl = false;
        std::vector<Extruder> tmp_extruder;
        for (const auto &i : m_vector_extr)
            // tmp 10000
            if (( m_koff * integrator_timestep * 10000 ) < dist(mt))
                tmp_extruder.push_back(*i);
	
	m_vector_extr.clear();
	for(auto &i : tmp_extruder)
		m_vector_extr.push_back(std::make_unique<Extruder>(i));

	if(m_extr.place_extruder(poly)){
            if(m_vector_extr.size() != 0) 
                for (auto &j : m_vector_extr)
                    if ((*j).extr_overlap(m_extr)) {
                        is_overl = true;
                        break;
                    }
        // tmp 10000
	    if (poly.get_poly_nmonomers() * m_kon * 0.032
			    * integrator_timestep * 10000 > dist(mt) and !(is_overl))
            	m_vector_extr.push_back(std::make_unique<Extruder>(m_extr));
    	}
    }
}

bool VectorExtruder::overlap_lr(Extruder &extr) {
    bool is_overl = false;
    for (auto &i : m_vector_extr) {
        // same extrusor
        if (extr == (*i))
            is_overl = false;

        else if (extr.get_l()-1 == (*i).get_r()) {
            is_overl = true;
            break;
        }
    }
    return is_overl;
}

bool VectorExtruder::overlap_ll(Extruder &extr) {
    bool is_overl = false;
    for (auto &i : m_vector_extr) {
        // same extrusor
        if (extr == (*i))
            is_overl = false;

        else if (extr.get_l()-1 == (*i).get_l()) {
            is_overl = true;
            break;
        }
    }
    return is_overl;
}
bool VectorExtruder::overlap_rl(Extruder &extr) {
    bool is_overl = false;
    for (auto &i : m_vector_extr) {
        // same extrusor
        if (extr == (*i))
            is_overl = false;

        else if (extr.get_r()+1 == (*i).get_l()) {
            is_overl = true;
            break;
        }
    }
    return is_overl;
}

bool VectorExtruder::overlap_rr(Extruder &extr) {
    bool is_overl = false;
    for (auto &i : m_vector_extr) {
        // same extrusor
        if (extr == (*i))
            is_overl = false;

        else if (extr.get_r()+1 == (*i).get_r()) {
            is_overl = true;
            break;
        }
    }
    return is_overl;
}
