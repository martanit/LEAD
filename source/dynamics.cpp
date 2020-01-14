#include "dynamics.h"
bool print_sys(Polymer &, VectorExtruder, std::string) ;

void Dynamics::run(bool rouse, bool soft_core, bool lennard_jones, bool compute_energy, std::string output) {
    for(unsigned long int i=0; i<m_dynamics_nstep; ++i) {
        if(compute_energy) m_poly_old = *m_poly;
        (*m_poly).reset_force();
        if(compute_energy) (*m_poly).reset_energy();

        Potential::set_new_polymer(*m_poly);

        this->box(compute_energy);
        if(rouse)
            this->harmonic_spring_f(compute_energy);
        else if(soft_core) {
            this->harmonic_spring_f(compute_energy);
            this->soft_core_f(i, compute_energy);
        }
        else if(lennard_jones) {
            this->harmonic_spring_f(compute_energy);
            this->lennard_jones_f(i, compute_energy);
        }

        m_poly = std::make_unique<Polymer>(Potential::get_poly());

        Integrator::set_new_polymer(*m_poly);

        this->langevin_overdamped();

        m_poly = std::make_unique<Polymer>(Integrator::get_poly());

        if(i%m_dynamics_print == 0) {
            print_xyz(*m_poly, output);
            if(compute_energy) std::cout << i*m_parm.get_timestep()/1E12 << " " << delta_h() << std::endl;
        }
    }
}

void Dynamics::run_extrusion(bool rouse, bool soft_core, bool lennard_jones, bool compute_energy, bool homogeneus_density, std::string output) {
    for (unsigned long int i = 0; i < m_dynamics_nstep; ++i) {
        if(compute_energy) m_poly_old = *m_poly;

        if (i % m_dynamics_print == 0) {
            if(homogeneus_density)
                m_vector_extr.update(*m_poly);
            else
                m_vector_extr.update_diff_density(*m_poly, i);
        }

        (*m_poly).reset_force();
        if(compute_energy) (*m_poly).reset_energy();

        Potential::set_new_polymer(*m_poly);
        Potential::set_new_extruder(m_vector_extr);

        this->box(compute_energy);
        this->extruder_spring_f(compute_energy);
        if(rouse)
            this->harmonic_spring_f(compute_energy);
        else if(soft_core) {
            this->harmonic_spring_f(compute_energy);
            this->soft_core_f(i, compute_energy);
        }
        else if(lennard_jones) {
            this->harmonic_spring_f(compute_energy);
            this->lennard_jones_f(i, compute_energy);
        }

        m_poly = std::make_unique<Polymer>(Potential::get_poly());
        m_vector_extr = Potential::get_extr();

        Integrator::set_new_polymer(*m_poly);
        Integrator::set_new_extruder(m_vector_extr);

        this->markov_chain();
        this->langevin_overdamped();

        m_poly = std::make_unique<Polymer>(Integrator::get_poly());
        m_vector_extr = Integrator::get_extr();

        if (i % m_dynamics_print == 0) {
            print_sys(output);
            if(compute_energy) std::cout << i*m_parm.get_timestep()/1E12 << " " << delta_h() << std::endl;
            int num_extr = 0;
            for (auto &i : m_vector_extr) {
                print_r(*m_poly, *i,
                        output + std::to_string(num_extr) + ".le");
                ++num_extr;
            }
        }

    }
}

void Dynamics::run_extrusion_field(bool rouse, bool soft_core, bool lennard_jones, bool compute_energy, std::string output) {
    for (unsigned long int i = 0; i < m_dynamics_nstep; ++i) {
        if(compute_energy) m_poly_old = *m_poly;
        if (i % m_dynamics_print == 0)
            m_vector_extr.update_field(*m_poly, m_interaction);
        (*m_poly).reset_force();
        if(compute_energy) (*m_poly).reset_energy();
        Potential::set_new_polymer(*m_poly);
        Potential::set_new_extruder(m_vector_extr);

        this->box(compute_energy);
        this->extruder_spring_f(compute_energy);
        if(rouse)
            this->harmonic_spring_f(compute_energy);
        else if(soft_core) {
            this->harmonic_spring_f(compute_energy);
            this->soft_core_f(i, compute_energy);
        }
        else if(lennard_jones) {
            this->harmonic_spring_f(compute_energy);
            this->lennard_jones_f(i, compute_energy);
        }

        m_poly = std::make_unique<Polymer>(Potential::get_poly());
        m_vector_extr = Potential::get_extr();

        Integrator::set_new_polymer(*m_poly);
        Integrator::set_new_extruder(m_vector_extr);
        Integrator::set_new_field(m_interaction);

        this->markov_chain();
        this->langevin_overdamped();
        this->extruders_diffusion();

        m_poly = std::make_unique<Polymer>(Integrator::get_poly());
        m_vector_extr = Integrator::get_extr();
        m_interaction = Integrator::get_field();

        if (i % m_dynamics_print == 0) {
            print_xyz(*m_poly, output);
            if(compute_energy) std::cout << i*m_parm.get_timestep()/1E12 << " " << delta_h() << std::endl;
            int num_extr = 0;
            for (auto &i : m_vector_extr) {
                print_r(*m_poly, *i,
                        output + std::to_string(num_extr) + ".le");
                ++num_extr;
            }
            print_field(m_interaction, "output/field.xyz");
        }
    }
}

double Dynamics::delta_h() {

    double first_term=0;
    double second_term=0;
    double delta_energy;

    for (unsigned int j=0; j< m_parm.get_nmonomers(); ++j) {
        first_term+=m_parm.get_timestep()/(4.*m_parm.get_gamma())*
                    ((std::pow((*m_poly).get_fx(j),2)+
                      std::pow((*m_poly).get_fy(j),2)+
                      std::pow((*m_poly).get_fz(j),2))-
                     (std::pow(m_poly_old.get_fx(j),2)+
                      std::pow(m_poly_old.get_fy(j),2)+
                      std::pow(m_poly_old.get_fz(j),2)));

        second_term+=(((*m_poly).get_x(j)-m_poly_old.get_x(j))*
                      ((*m_poly).get_fx(j)+m_poly_old.get_fx(j))/2.+
                      ((*m_poly).get_y(j)-m_poly_old.get_y(j))*
                      ((*m_poly).get_fy(j)+m_poly_old.get_fy(j))/2.+
                      ((*m_poly).get_z(j)-m_poly_old.get_z(j))*
                      ((*m_poly).get_fz(j)+m_poly_old.get_fz(j))/2.);
    }

    delta_energy=(*m_poly).get_energy()-m_poly_old.get_energy();

    return first_term+second_term+delta_energy;
}

bool Dynamics::print_sys(std::string out_xyz) {
    // open stream to write xyz file
    std::ofstream output;
    output.open(out_xyz, std::fstream::app);

    // return error if read file fail
    if (output.fail()) {
        throw "ERROR: Impossible to write xyz file " + out_xyz;
        return 1;
    }


    std::vector<bool> v((*m_poly).get_poly_nmonomers());
    if(m_vector_extr.vextr_size() !=0) {
        for (int i = 0; i < (*m_poly).get_poly_nmonomers(); i++) {
            for(auto &a : m_vector_extr) {
                if((*a).get_l() == i) {
                    v[i] = 1;
                    break;
                }
                else if((*a).get_r() == i) {
                    v[i] = 1;
                    break;
                }
                else v[i] = 0;
            }
        }
    }
    else std::fill(v.begin(), v.end(), 0);

    output << (*m_poly).get_poly_nmonomers() << std::endl << std::endl;
    for (int i = 0; i < (*m_poly).get_poly_nmonomers(); i++) {
        if(v[i] == 1)
            output << "\tLe" << "\t\t\t"
                   << (*m_poly).get_x(i) << "\t\t\t"
                   << (*m_poly).get_y(i) << "\t\t\t"
                   << (*m_poly).get_z(i) << std::endl;

        else
            output << "\tAu" << "\t\t\t"
                   << (*m_poly).get_x(i) << "\t\t\t"
                   << (*m_poly).get_y(i) << "\t\t\t"
                   << (*m_poly).get_z(i) << std::endl;
    }
    output.close();
    return 0;
}

