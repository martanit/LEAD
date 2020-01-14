#include "potential.h"

void Potential::lennard_jones_f(int step,
                                bool compute_energy) {
    if (step % 10000 == 0 or step == 0) {
        for (unsigned int i = 0; i < m_poly.get_poly_nmonomers(); ++i) {
            sphere[i].clear();
            for (unsigned int j = i + 1; j < m_poly.get_poly_nmonomers(); ++j) {
                x = (m_poly.get_x(i) - m_poly.get_x(j));
                y = (m_poly.get_y(i) - m_poly.get_y(j));
                z = (m_poly.get_z(i) - m_poly.get_z(j));

                dr = x * x + y * y + z * z;

                // calculate points that are into m_pot_rcut^2
                if (dr < m_pot_rcut * m_pot_rcut) {
                    sphere[i].push_back(j);
                }
            }
        }
    }
    for (unsigned int i = 0; i < m_poly.get_poly_nmonomers(); ++i) {
        for (auto &&k : sphere[i]) {
            x = (m_poly.get_x(i) - m_poly.get_x(k));
            y = (m_poly.get_y(i) - m_poly.get_y(k));
            z = (m_poly.get_z(i) - m_poly.get_z(k));

            dr = x * x + y * y + z * z;
            if (dr < m_pot_rcut * m_pot_rcut) {
                f_x = x * m_pot_epsilon * 12.0 * m_pot_rmin_12 / std::pow(dr, 7);
                f_y = y * m_pot_epsilon * 12.0 * m_pot_rmin_12 / std::pow(dr, 7);
                f_z = z * m_pot_epsilon * 12.0 * m_pot_rmin_12 / std::pow(dr, 7);

                f_x -= x * m_pot_epsilon * 12.0 * m_pot_rmin_6 / std::pow(dr, 4);
                f_y -= y * m_pot_epsilon * 12.0 * m_pot_rmin_6 / std::pow(dr, 4);
                f_z -= z * m_pot_epsilon * 12.0 * m_pot_rmin_6 / std::pow(dr, 4);

                m_poly.add_force(i, f_x, f_y, f_z);
                m_poly.add_force(k, -f_x, -f_y, -f_z);

                if (compute_energy)
                    m_poly.add_energy(-2. * m_pot_epsilon * m_pot_rmin_6 /
                                      std::pow(dr, 3) + m_pot_epsilon * m_pot_rmin_12 / std::pow(dr, 6));
            }
        }
    }
}

void Potential::soft_core_f(int step,
                            bool compute_energy) {
    if (step % 10000 == 0 or step == 0) {
        for (unsigned int i = 0; i < m_poly.get_poly_nmonomers(); ++i) {
            sphere[i].clear();
            for (unsigned int j = i + 1; j < m_poly.get_poly_nmonomers(); ++j) {
                x = (m_poly.get_x(i) - m_poly.get_x(j));
                y = (m_poly.get_y(i) - m_poly.get_y(j));
                z = (m_poly.get_z(i) - m_poly.get_z(j));

                dr = x * x + y * y + z * z;

                // calculate points that are into m_pot_rcut^2
                if (dr < m_pot_rcut * m_pot_rcut) {
                    sphere[i].push_back(j);
                }
            }
        }
    }

    for (unsigned int i = 0; i < m_poly.get_poly_nmonomers(); ++i) {
        for (auto &&k : sphere[i]) {

            x = (m_poly.get_x(i) - m_poly.get_x(k));
            y = (m_poly.get_y(i) - m_poly.get_y(k));
            z = (m_poly.get_z(i) - m_poly.get_z(k));

            dr = x * x + y * y + z * z;
            if (dr < m_pot_rcut * m_pot_rcut) {
                f_x = x * m_pot_epsilon * 12.0 * m_pot_rmin_12 / std::pow(dr, 7);
                f_y = y * m_pot_epsilon * 12.0 * m_pot_rmin_12 / std::pow(dr, 7);
                f_z = z * m_pot_epsilon * 12.0 * m_pot_rmin_12 / std::pow(dr, 7);

                m_poly.add_force(i, f_x, f_y, f_z);
                m_poly.add_force(k, -f_x, -f_y, -f_z);

                if (compute_energy)
                    m_poly.add_energy(m_pot_epsilon * m_pot_rmin_12 / std::pow(dr, 6));
            }
        }
    }
}

void Potential::harmonic_spring_f(bool compute_energy) {
    for (int i = 0; i < m_poly.get_poly_nmonomers() - 1; ++i) {
        spring_x = -k * (m_poly.dist(i, i + 1) - m_pot_rmin) *
                   (m_poly.get_x(i) - m_poly.get_x(i + 1)) / m_poly.dist(i, i + 1);
        spring_y = -k * (m_poly.dist(i, i + 1) - m_pot_rmin) *
                   (m_poly.get_y(i) - m_poly.get_y(i + 1)) / m_poly.dist(i, i + 1);
        spring_z = -k * (m_poly.dist(i, i + 1) - m_pot_rmin) *
                   (m_poly.get_z(i) - m_poly.get_z(i + 1)) / m_poly.dist(i, i + 1);

        m_poly.add_force(i, spring_x, spring_y, spring_z);
        m_poly.add_force(i + 1, -spring_x, -spring_y, -spring_z);
        if(compute_energy) m_poly.add_energy(k*(m_poly.dist(i,i+1)-m_pot_rmin)*
                                                 (m_poly.dist(i,i+1)-m_pot_rmin)/2.);
    }
}

void Potential::extruder_spring_f(bool compute_energy) {
    // iterate over each extruder
    for (auto &i : m_vector_extr) {
        spring_x = -k_extr *
                   (m_poly.dist((*i).get_l(), (*i).get_r()) - extr_length) *
                   (m_poly.get_x((*i).get_l()) - m_poly.get_x((*i).get_r())) /
                   m_poly.dist((*i).get_l(), (*i).get_r());
        spring_y = -k_extr *
                   (m_poly.dist((*i).get_l(), (*i).get_r()) - extr_length) *
                   (m_poly.get_y((*i).get_l()) - m_poly.get_y((*i).get_r())) /
                   m_poly.dist((*i).get_l(), (*i).get_r());
        spring_z = -k_extr *
                   (m_poly.dist((*i).get_l(), (*i).get_r()) - extr_length) *
                   (m_poly.get_z((*i).get_l()) - m_poly.get_z((*i).get_r())) /
                   m_poly.dist((*i).get_l(), (*i).get_r());

        m_poly.add_force((*i).get_l(), spring_x, spring_y, spring_z);
        m_poly.add_force((*i).get_r(), -spring_x, -spring_y, -spring_z);

        if(compute_energy) m_poly.add_energy(k*(m_poly.dist((*i).get_l(),(*i).get_r())-extr_length)*
                                                 (m_poly.dist((*i).get_l(),(*i).get_r())-extr_length)/2.);
    }
}

void Potential::box(bool compute_energy) {
    for (unsigned int i = 0; i < m_poly.get_poly_nmonomers(); ++i) {
        if( m_poly.get_x(i) < -m_pot_box_length/2. and
                m_poly.get_y(i) < -m_pot_box_length/2. and
                m_poly.get_z(i) < -m_pot_box_length/2. ) {
            f_x = -m_pot_kside *(m_poly.get_x(i)+m_pot_box_length/2.);
            f_y = -m_pot_kside *(m_poly.get_y(i)+m_pot_box_length/2.);
            f_z = -m_pot_kside *(m_poly.get_z(i)+m_pot_box_length/2.);

            m_poly.add_force(i, f_x, f_y, f_z);

            if(compute_energy) m_poly.add_energy(m_pot_kside/2.*
                                                     (std::pow(m_poly.get_x(i)+m_pot_box_length/2., 2)+
                                                      std::pow(m_poly.get_y(i)+m_pot_box_length/2., 2)+
                                                      std::pow(m_poly.get_z(i)+m_pot_box_length/2., 2)));
        }
        else if( m_poly.get_x(i) > m_pot_box_length/2. and
                 m_poly.get_y(i) > m_pot_box_length/2. and
                 m_poly.get_z(i) > m_pot_box_length/2. ) {
            f_x = -m_pot_kside *(m_poly.get_x(i)-m_pot_box_length/2.);
            f_y = -m_pot_kside *(m_poly.get_y(i)-m_pot_box_length/2.);
            f_z = -m_pot_kside *(m_poly.get_z(i)-m_pot_box_length/2.);

            m_poly.add_force(i, f_x, f_y, f_z);

            if(compute_energy) m_poly.add_energy(m_pot_kside/2.*
                                                     (std::pow(m_poly.get_x(i)+m_pot_box_length/2., 2)+
                                                      std::pow(m_poly.get_y(i)+m_pot_box_length/2., 2)+
                                                      std::pow(m_poly.get_z(i)+m_pot_box_length/2., 2)));
        }
    }
}
