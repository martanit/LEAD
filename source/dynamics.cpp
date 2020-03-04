#include "dynamics.h"
bool print_sys(Polymer &, VectorExtruder, std::string) ;

void Dynamics::run(bool rouse, bool soft_core, bool lennard_jones, bool compute_energy, std::string output) {
    for(unsigned long int i=0; i<m_dynamics_nstep; ++i) {

        if(compute_energy) m_poly_old = Integrator::m_poly;

	Integrator::m_poly.reset_force();
        if(compute_energy) Integrator::m_poly.reset_energy();

        Potential::m_poly = Integrator::m_poly;

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
	    this->well(compute_energy);
        }

        Integrator::m_poly = Potential::m_poly;

        this->langevin_overdamped();

        if(i%m_dynamics_print == 0) {
            print_xyz(Integrator::m_poly, output);
            if(compute_energy) std::cout << i*m_parm.get_timestep()/1E12 << " " << delta_h() << std::endl;
        }
    }
}

void Dynamics::run_extrusion(bool rouse, bool soft_core, bool lennard_jones, bool compute_energy, std::string output) {
    std::ofstream nextr_out;
    nextr_out.open(output + ".le", std::fstream::app);
    if (nextr_out.fail()) 
            throw "ERROR: Impossible to write number of extruder to "
                                                     + output + ".le";

    for (unsigned long int i = 0; i < m_dynamics_nstep; ++i) {
        
    	if(compute_energy) m_poly_old = Integrator::m_poly;
        if (i % m_dynamics_print == 0) 
    		Integrator::m_vector_extr.update(Integrator::m_poly);
    
    	Integrator::m_poly.reset_force();
        if(compute_energy) Integrator::m_poly.reset_energy();
    
        Potential::m_poly = Integrator::m_poly;
        if (i % m_dynamics_print == 0) 
            Potential::m_vector_extr = Integrator::m_vector_extr;
    
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
    
        Integrator::m_poly = Potential::m_poly;
        if (i % m_dynamics_print == 0) 
    	    Integrator::m_vector_extr = Potential::m_vector_extr;
    
        if (i % m_dynamics_print == 0) 
            this->markov_chain();
        this->langevin_overdamped();
    
        if (i % m_dynamics_print == 0) {
            print_sys(output);
            if(compute_energy) std::cout << i*m_parm.get_timestep()/1E12 
                                         << " " << delta_h() << std::endl;
            int num_extr = 0;
            for (auto &i : Integrator::m_vector_extr) {
           //     print_r(*i,
           //             output + std::to_string(num_extr) + ".le");
                ++num_extr;
            }
            nextr_out << i << " " << num_extr << std::endl;
        }
    }
}


double Dynamics::delta_h() {

    double first_term=0;
    double second_term=0;
    double delta_energy;

    for (unsigned int j=0; j< Potential::m_poly.get_poly_nmonomers(); ++j) {
        first_term+=m_parm.get_timestep()/(4.*m_parm.get_gamma())*
                    ((std::pow(Integrator::m_poly.get_fx(j),2)+
                      std::pow(Integrator::m_poly.get_fy(j),2)+
                      std::pow(Integrator::m_poly.get_fz(j),2))-
                     (std::pow(m_poly_old.get_fx(j),2)+
                      std::pow(m_poly_old.get_fy(j),2)+
                      std::pow(m_poly_old.get_fz(j),2)));

        second_term+=((Integrator::m_poly.get_x(j)-m_poly_old.get_x(j))*
                      (Integrator::m_poly.get_fx(j)+m_poly_old.get_fx(j))/2.+
                      (Integrator::m_poly.get_y(j)-m_poly_old.get_y(j))*
                      (Integrator::m_poly.get_fy(j)+m_poly_old.get_fy(j))/2.+
                      (Integrator::m_poly.get_z(j)-m_poly_old.get_z(j))*
                      (Integrator::m_poly.get_fz(j)+m_poly_old.get_fz(j))/2.);
    }

    delta_energy=Integrator::m_poly.get_energy()-m_poly_old.get_energy();

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


    std::vector<bool> v(Integrator::m_poly.get_poly_nmonomers());
    if(Integrator::m_vector_extr.vextr_size() !=0) {
        for (int i = 0; i < Integrator::m_poly.get_poly_nmonomers(); i++) {
            for(auto &a : Integrator::m_vector_extr) {
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

    output << Integrator::m_poly.get_poly_nmonomers() << std::endl << std::endl;
    for (int i = 0; i < Integrator::m_poly.get_poly_nmonomers(); i++) {
        if(v[i] == 1)
            output << "\tLe" << "\t\t\t"
                   << Integrator::m_poly.get_x(i) << "\t\t\t"
                   << Integrator::m_poly.get_y(i) << "\t\t\t"
                   << Integrator::m_poly.get_z(i) << std::endl;

        else
            output << "\tAu" << "\t\t\t"
                   << Integrator::m_poly.get_x(i) << "\t\t\t"
                   << Integrator::m_poly.get_y(i) << "\t\t\t"
                   << Integrator::m_poly.get_z(i) << std::endl;
    }
    output.close();
    return 0;
}

