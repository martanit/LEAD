#include "vector_extruder.h"

void VectorExtruder::first_fill(Polymer &poly) {
    bool is_overl = false;

    std::vector<Extruder> tmp_extruder;
    for (int i = 0; i < m_n_max_extr; ++i) {
        m_extr.place_extruder(poly);
        if( tmp_extruder.size() != 0) {
            for (auto &j : tmp_extruder)
                if (j.extr_overlap(m_extr)) {
                    is_overl = true;
                    break;
                }
        }
        if ((m_kon * integrator_timestep /
                poly.get_poly_nmonomers()) > dist(mt) and
                !(is_overl))
            tmp_extruder.push_back(m_extr);
    }

    m_vector_extr.clear();
    for (const auto &i : tmp_extruder)
        m_vector_extr.push_back(std::make_unique<Extruder>(i));
}

void VectorExtruder::update(Polymer &poly) {

    if(m_vector_extr.size() == 0) this->first_fill(poly);

    else {
        bool is_overl = false;
        std::vector<Extruder> tmp_extruder;
        for (const auto &i : m_vector_extr)
            // fill tmp_extruder only with extruder
            // that are not undbind
            if ((m_vector_extr.size() * m_koff * integrator_timestep) < dist(mt))
                tmp_extruder.push_back(*i);

        for (int i = 0; i < m_n_max_extr - tmp_extruder.size(); ++i) {
            m_extr.place_extruder(poly);
            for (auto &j : tmp_extruder)
                if ((j).extr_overlap(m_extr)) {
                    is_overl = true;
                    break;
                }
            if ((m_kon * integrator_timestep / poly.get_poly_nmonomers()) >
                    dist(mt) and
                    !(is_overl))
                tmp_extruder.push_back(m_extr);
        }
        m_vector_extr.clear();
        for (const auto &i : tmp_extruder)
            m_vector_extr.push_back(std::make_unique<Extruder>(i));
    }
}

void VectorExtruder::update_diff_density(Polymer &poly, unsigned long int t) {

    if(m_vector_extr.size() == 0) this->first_fill(poly);

    else {
        bool is_overl = false;
        std::vector<Extruder> tmp_extruder;
        for (const auto &i : m_vector_extr) {
            // fill tmp_extruder only with extruder
            // that are not undbind, save unbind extruders
            if ((m_vector_extr.size() * m_koff * integrator_timestep) < dist(mt))
                tmp_extruder.push_back(*i);
	    else 
		 m_unloaded_extr.push_back({*i, t});
	}
	
	// Move unloaded extruders from unbinding vector to background
	m_unloaded_extr.erase(std::remove_if(m_unloaded_extr.begin(), m_unloaded_extr.end(),
		[this,t](IdxTime x) {
			return x.t < t-dt; 
			}), m_unloaded_extr.end());

        for (int i = 0; i < m_n_max_extr - tmp_extruder.size(); ++i) {
	    m_extr.place_extruder(poly);
	    r = m_extr.xyz_position(poly);

            for (auto &j : tmp_extruder)
                if ((j).extr_overlap(m_extr)) {
                    is_overl = true;
                    break;
                }
            if ((m_kon/m_norm * this->density(r, t, poly) * integrator_timestep / poly.get_poly_nmonomers()) >
                    dist(mt) and
                    !(is_overl))
                tmp_extruder.push_back(m_extr);
        }
        m_vector_extr.clear();
        for (const auto &i : tmp_extruder)
            m_vector_extr.push_back(std::make_unique<Extruder>(i));
    }
}

double VectorExtruder::density(Position r, unsigned long int t, Polymer &poly) {
	double sum=0;
	// Number of unloaded extruder is total 
	// extruder minus bounded extrduer
	Nu = N0-m_vector_extr.size(); 
	for (auto &i : m_unloaded_extr) {
	    ri = i.extr.xyz_position(poly);
	    sum = sum + 1./(V*std::sqrt((2*M_PI*D*(t-i.t))))*std::exp(-1*(std::pow((r.x-ri.x),2)+
				     std::pow((r.y-ri.y),2)+
				     std::pow((r.z-ri.z),2))/(
				     2.*(t-i.t)*D));
	}
	return (Nu-m_unloaded_extr.size())/V+sum;
}

bool VectorExtruder::overlap_lr(Extruder &extr) {
    bool is_overl = false;
    for (auto &i : m_vector_extr) {
        // same extrusor
        if (extr == (*i))
            is_overl = false;

        else if (extr.get_l() == (*i).get_r()) {
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

        else if (extr.get_l() == (*i).get_l()) {
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

        else if (extr.get_r() == (*i).get_l()) {
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

        else if (extr.get_r() == (*i).get_r()) {
            is_overl = true;
            break;
        }
    }
    return is_overl;
}
