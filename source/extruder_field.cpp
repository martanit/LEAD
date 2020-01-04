#include "extruder_field.h"

ExtruderField::ExtruderField(Parameters parm)
      :m_field_step(parm.get_field_step()),
      m_k_diff(parm.get_k_diff()),
      uniform01(0., 1.),
      uniform05(0, 5),
      m_extruder_c(m_field_length, std::vector<std::vector<double>>(m_field_length, std::vector<double>(m_field_length))){
	this->first_extruders_field();
}

// assuming uniform distribution for extruder concentration
void ExtruderField::first_extruders_field() {
    for(auto &i : m_extruder_c)
        for(auto &j : i)
            for(auto &k : j)
                k = uniform01(mt);
}

void ExtruderField::add_delta_c(int i, int j, int k) {
	m_extruder_c[i][j][k] += m_delta_c;
}

