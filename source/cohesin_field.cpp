#include "cohesin_field.h"

CohesinField::CohesinField(Parameters parm)
    : m_poly_sphere(parm.get_psphere()),
      m_field_step(parm.get_field_step()),
      m_k_diff(parm.get_k_diff()),
      uniform01(0., 1.),
      uniform05(0, 5),
      m_cohesin_c(m_poly_sphere, std::vector<std::vector<double>>(m_poly_sphere, std::vector<double>(m_poly_sphere))),
      m_cohesin_c_new(m_poly_sphere, std::vector<std::vector<double>>(m_poly_sphere, std::vector<double>(m_poly_sphere))){
	this->set_cohesin();
}

CohesinField::~CohesinField() {  };


// assuming uniform distribution for cohesin concentration
void CohesinField::set_cohesin() {
    for(auto &i : m_cohesin_c)
        for(auto &j : i)
            for(auto &k : j)
                k = uniform01(mt);
}


