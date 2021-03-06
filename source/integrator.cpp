#include "integrator.h"

void Integrator::langevin_overdamped() {
    for (unsigned int i = 0; i < m_poly.get_poly_nmonomers(); ++i) {
        m_poly.set_x(m_poly.get_x(i) +
                        (m_poly.get_fx(i) + stoch_term * gauss_term(mt)) *
                        m_integrator_timestep / m_integrator_gamma,
                        i);

        m_poly.set_y(m_poly.get_y(i) +
                        (m_poly.get_fy(i) + stoch_term * gauss_term(mt)) *
                        m_integrator_timestep / m_integrator_gamma,
                        i);

        m_poly.set_z(m_poly.get_z(i) +
                        (m_poly.get_fz(i) + stoch_term * gauss_term(mt)) *
                        m_integrator_timestep / m_integrator_gamma,
                        i);
    }
}

void Integrator::markov_chain() {
    for (auto &&i : m_vector_extr) {
        // right side go foward to right
        if ((*i).get_r() != (m_poly.get_poly_nmonomers() - 1))
            if((((*i).get_ctcf()[(*i).get_r() + 1] >= 0 and 
                 (*i).get_rate_fwr() * m_integrator_timestep * m_integrator_deltastep > 
                    transition_prob(mt)) or 
                ((*i).get_ctcf()[(*i).get_r() + 1] < 0 and
                 (*i).get_rate_fwr() * (*i).get_permeability() /
                    std::abs((*i).get_ctcf()[(*i).get_r() + 1]) * 
                    m_integrator_timestep * m_integrator_deltastep > transition_prob(mt))) and
                    !(m_vector_extr.overlap_rl(*i)) and
		    !(m_vector_extr.overlap_rr(*i)))
                (*i).set_r((*i).get_r() + 1);
        // right side go backward to left
        if ((((*i).get_ctcf()[(*i).get_r() - 1] <= 0 and
              (*i).get_rate_bwr() * m_integrator_timestep * m_integrator_deltastep > 
                 transition_prob(mt)) or
             ((*i).get_ctcf()[(*i).get_r() - 1] > 0 and
                 (*i).get_rate_bwr() * (*i).get_permeability() /
                 std::abs((*i).get_ctcf()[(*i).get_r() - 1]) *
                 m_integrator_timestep * m_integrator_deltastep > transition_prob(mt))) and
                !(m_vector_extr.overlap_rr(*i)))
            (*i).set_r((*i).get_r() - 1);

        // left side go forward to left
        if ((*i).get_l() != 0)
            if ((((*i).get_ctcf()[(*i).get_l() - 1] <= 0 and
                  (*i).get_rate_fwl() * m_integrator_timestep * m_integrator_deltastep >
                    transition_prob(mt)) or
                 ((*i).get_ctcf()[(*i).get_l() - 1] > 0 and
                  (*i).get_rate_fwl() * (*i).get_permeability() /
                     std::abs((*i).get_ctcf()[(*i).get_l() - 1]) * 
                     m_integrator_timestep * m_integrator_deltastep > transition_prob(mt))) and
                    !(m_vector_extr.overlap_lr(*i)) and
		    !(m_vector_extr.overlap_ll(*i)))
                (*i).set_l((*i).get_l() - 1);
        // left side go backward to right
        if ((((*i).get_ctcf()[(*i).get_l() + 1] >= 0 and
              (*i).get_rate_bwl() * m_integrator_timestep * m_integrator_deltastep > 
                 transition_prob(mt)) or
             ((*i).get_ctcf()[(*i).get_l() + 1] < 0 and
              (*i).get_rate_bwl() * (*i).get_permeability() /
                 std::abs((*i).get_ctcf()[(*i).get_l() + 1]) *
                 m_integrator_timestep * m_integrator_deltastep >  transition_prob(mt))) and
                !(m_vector_extr.overlap_ll(*i)))
            (*i).set_l((*i).get_l() + 1);
    }
}
