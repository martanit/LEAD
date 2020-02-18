#include "vector_extruder.h"

void VectorExtruder::first_fill(Polymer &poly) {
    bool is_overl = false;
    bool is_set;
    is_set = m_extr.place_extruder(poly);
    if(is_set){
        if( m_vector_extr.size() != 0) {
            for (auto &j : m_vector_extr)
                if ((*j).extr_overlap(m_extr)) {
                    is_overl = true;
                    break;
                }
        }
        //tmp 10000
        if ((10000*poly.get_poly_nmonomers() * m_kon * 0.032 *
             integrator_timestep) > dist(mt) and !(is_overl))
            m_vector_extr.push_back(std::make_unique<Extruder>(m_extr));
    }
}

void VectorExtruder::first_fill_field(Polymer &poly, FieldAction cohes_field_int) {
    bool is_overl = false;
    int extr_per_cell;
    std::vector<Extruder> tmp_extruder;
    std::vector<Extruder> tmp_extruder_cell;

    cohes_field_int.interaction();

    for(auto &a : cohes_field_int.get_contact_cell()) {
        if(m_extr.can_place_extr(poly,cohes_field_int.monomer_min(a),
                                 cohes_field_int.monomer_max(a))) {
            tmp_extruder_cell.clear();
            for( int i = 0; i< cohes_field_int.get_c(a.i, a.j, a.k); ++i) {
                m_extr.place_extruder_cell(poly, cohes_field_int.monomer_min(a), cohes_field_int.monomer_max(a));
                if(tmp_extruder_cell.size()!=0)
                    for (auto &j : tmp_extruder_cell)
                        if (j.extr_overlap(m_extr)) {
                            is_overl = true;
                            break;
                        }
                if (( m_kon*integrator_timestep/cohes_field_int.get_init_c()) >
                        dist(mt) and !(is_overl)) {
                    tmp_extruder_cell.push_back(m_extr);
                    cohes_field_int.sub_delta_c(a.i,a.j,a.k);
                }
            }
            for(auto &e : tmp_extruder_cell)
                tmp_extruder.push_back(e);
        }
    }

    m_vector_extr.clear();
    for (const auto &i : tmp_extruder)
        m_vector_extr.push_back(std::make_unique<Extruder>(i));
}

void VectorExtruder::update(Polymer &poly) {

    if(m_vector_extr.size() == 0) this->first_fill(poly);

    else {
        bool is_overl = false;
        bool is_set;
        std::vector<Extruder> tmp_extruder;
        for (const auto &i : m_vector_extr) 
            //tmp 10000
            if(10000*m_koff * integrator_timestep < dist(mt))
                tmp_extruder.push_back(*i);
        
        m_vector_extr.clear();
        for (const auto &i : tmp_extruder)
            m_vector_extr.push_back(std::make_unique<Extruder>(i));
        
        is_set=m_extr.place_extruder(poly);
        
        if(is_set){
            for (auto &j : m_vector_extr)
                if ((*j).extr_overlap(m_extr)) {
                    is_overl = true;
                    break;
                }
            //tmp 10000
            if ((10000*poly.get_poly_nmonomers() * m_kon * 0.032 *
                 integrator_timestep) > dist(mt) and !(is_overl))
                m_vector_extr.push_back(std::make_unique<Extruder>(m_extr));
        }
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
            if (( m_koff * integrator_timestep) < dist(mt))
                tmp_extruder.push_back(*i);
            else
                m_unloaded_extr.push_back({*i, t});
        }

        // Move unloaded extruders from unbinding vector to background
        m_unloaded_extr.erase(std::remove_if(m_unloaded_extr.begin(), m_unloaded_extr.end(),
        [this,t](IdxTime x) {
            return x.t < t-m_dt;
        }), m_unloaded_extr.end());

        int tmp_size = tmp_extruder.size();
        for (int i = 0; i < m_n_max_extr - tmp_size; ++i) {
            m_extr.place_extruder(poly);
            r = m_extr.xyz_position(poly);

            for (auto &j : tmp_extruder)
                if ((j).extr_overlap(m_extr)) {
                    is_overl = true;
                    break;
                }
            if ((m_kon * this->density(r, t, poly)/rho0_tot * integrator_timestep / m_n_max_extr-tmp_size) >
                    dist(mt) and
                    !(is_overl))
                tmp_extruder.push_back(m_extr);
        }

        m_vector_extr.clear();
        for (const auto &i : tmp_extruder)
            m_vector_extr.push_back(std::make_unique<Extruder>(i));
    }
}

void VectorExtruder::update_field(Polymer &poly, FieldAction cohes_field_int) {

    if(m_vector_extr.size() == 0) this->first_fill_field(poly, cohes_field_int);
    else {
        bool is_overl = false;
        int extr_per_cell;
        std::vector<Extruder> tmp_extruder;
        std::vector<Extruder> tmp_extruder_cell;

        cohes_field_int.interaction();
        
        for (const auto &i : m_vector_extr)
            // fill tmp_extruder only with extruder
            // that are not undbind
            if (( m_koff * integrator_timestep) < dist(mt))
                tmp_extruder.push_back(*i);
            else {
                // Convention, we take the right monomer whose the extruder is bind as reference
                Cell a = cohes_field_int.monomer_cell(poly.get_x((*i).get_r()), poly.get_y((*i).get_r()), poly.get_z((*i).get_r()));
                cohes_field_int.add_delta_c(a.i,a.j,a.k);
            }
        
        for(auto &a : cohes_field_int.get_contact_cell()) {
            if(m_extr.can_place_extr(poly,cohes_field_int.monomer_min(a), cohes_field_int.monomer_max(a))) {
                tmp_extruder_cell.clear();
                for( int i = 0; i< cohes_field_int.get_c(a.i, a.j, a.k); ++i) {
                    m_extr.place_extruder_cell(poly, cohes_field_int.monomer_min(a), cohes_field_int.monomer_max(a));
                    if(tmp_extruder_cell.size()!=0)
                        for (auto &j : tmp_extruder_cell)
                            if (j.extr_overlap(m_extr)) {
                                is_overl = true;
                                break;
                            }
                    if (( m_kon*integrator_timestep/cohes_field_int.get_init_c()) >
                            dist(mt) and !(is_overl)){
                        tmp_extruder_cell.push_back(m_extr);
                        cohes_field_int.sub_delta_c(a.i,a.j,a.k);
                    }
                }
                for(auto &e : tmp_extruder_cell)
                    tmp_extruder.push_back(e);
            }
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
    Nb_eq = m_vector_extr.size();
    for (auto &i : m_unloaded_extr) {
        ri = i.extr.xyz_position(poly);
        sum = sum + 1./(poly.get_Veq()*std::sqrt((2*M_PI*D*(t-i.t))))*
              std::exp(-1*(std::pow((r.x-ri.x),2)+
                           std::pow((r.y-ri.y),2)+
                           std::pow((r.z-ri.z),2))/(
                           2.*(t-i.t)*D));
    }
    return rho0_tot-(Nb_eq+m_unloaded_extr.size())
           /poly.get_Veq()+sum;
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
