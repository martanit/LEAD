#include "integrator.h"

// integrator

void Integrator::langevin_overdamped() {
    for (unsigned int i = 0; i < (*m_poly).get_poly_nmonomers(); ++i) {
        (*m_poly).set_x((*m_poly).get_x(i) +
                        ((*m_poly).get_fx(i) + stoch_term * gauss_term(mt)) *
                        m_integrator_timestep / m_integrator_gamma,
                        i);

        (*m_poly).set_y((*m_poly).get_y(i) +
                        ((*m_poly).get_fy(i) + stoch_term * gauss_term(mt)) *
                        m_integrator_timestep / m_integrator_gamma,
                        i);

        (*m_poly).set_z((*m_poly).get_z(i) +
                        ((*m_poly).get_fz(i) + stoch_term * gauss_term(mt)) *
                        m_integrator_timestep / m_integrator_gamma,
                        i);
    }
}

void Integrator::markov_chain() {
    for (auto &&i : m_vector_extr) {
        // right side go foward to right
        if ((*i).get_r() != ((*m_poly).get_poly_nmonomers() - 1))
            if (((*i).get_rate_fwr() * m_integrator_timestep >
                    transition_prob(mt)) and
                    ((*i).get_ctcf()[(*i).get_r()] <= 0 or
                    (*i).get_rate_fwr() * (*i).get_permeability() /
                    (*i).get_ctcf()[(*i).get_r()] * m_integrator_timestep > 
                    transition_prob(mt)) and
                    !(m_vector_extr.overlap_rl(*i)))
                (*i).set_r((*i).get_r() + 1);
        // right side go backward to left
        if (((*i).get_rate_bwr() * m_integrator_timestep > transition_prob(mt)) and
                ((*i).get_ctcf()[(*i).get_r()] >= 0 or
                 (*i).get_rate_bwr() * (*i).get_permeability() /
                 (*i).get_ctcf()[(*i).get_r()] * m_integrator_timestep > 
                 transition_prob(mt)) and
                !(m_vector_extr.overlap_rr(*i)))
            (*i).set_r((*i).get_r() - 1);

        // left side go forward to left
        if ((*i).get_l() != 0)
            if (((*i).get_rate_fwl() * m_integrator_timestep >
                    transition_prob(mt)) and
                    ((*i).get_ctcf()[(*i).get_l()] >= 0 or
                     (*i).get_rate_fwl() * (*i).get_permeability() /
                     (*i).get_ctcf()[(*i).get_l()] * m_integrator_timestep > 
                     transition_prob(mt)) and
                    !(m_vector_extr.overlap_lr(*i)))
                (*i).set_l((*i).get_l() - 1);
        // left side go backward to right
        if (((*i).get_rate_bwl() * m_integrator_timestep > transition_prob(mt)) and
                ((*i).get_ctcf()[(*i).get_l()] <= 0 or
                 (*i).get_rate_bwl() * (*i).get_permeability() /
                 (*i).get_ctcf()[(*i).get_l()] * m_integrator_timestep > 
                 transition_prob(mt)) and
                !(m_vector_extr.overlap_ll(*i)))
            (*i).set_l((*i).get_l() + 1);
    }
}

void Integrator::extruders_diffusion() {
    bool move;
    for(int i = 0; i< m_field.get_field_length(); ++i)
        for(int j = 0; j < m_field.get_field_length(); ++j)
            for (int k = 0; k < m_field.get_field_length(); ++k) {
                if(m_field.get_c(i,j,k)*m_field.get_k_diff()*
                        m_integrator_timestep/m_field.get_field_step() > transition_prob(mt)
                        and m_field.get_c(i,j,k)>=m_field.get_delta_c()) {

                    move = false;
                    while(move != true) {
                        direction = uniform05(mt);
                        if(direction == 0 and i != 0) {
                            move = true;
                            new_field.add_delta_c(i-1,j,k);
                            new_field.sub_delta_c(i,j,k);
                        }
                        else if(direction == 1 and i != m_field.get_field_length()-1) {
                            move = true;
                            new_field.add_delta_c(i+1,j,k);
                            new_field.sub_delta_c(i,j,k);
                        }

                        else if(direction == 2 and j != 0) {
                            move = true;
                            new_field.add_delta_c(i,j-1,k);
                            new_field.sub_delta_c(i,j,k);
                        }
                        else if(direction == 3 and j != m_field.get_field_length()-1) {
                            move = true;
                            new_field.add_delta_c(i,j+1,k);
                            new_field.sub_delta_c(i,j,k);
                        }
                        else if(direction == 4 and k != 0) {
                            move = true;
                            new_field.add_delta_c(i,j,k-1);
                            new_field.sub_delta_c(i,j,k);
                        }
                        else if(direction == 5 and k != m_field.get_field_length()-1) {
                            move = true;
                            new_field.add_delta_c(i,j,k+1);
                            new_field.sub_delta_c(i,j,k);
                        }
                    }
                }
            }
    m_field = new_field;
}
