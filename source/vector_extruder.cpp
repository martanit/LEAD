#include "vector_extruder.h"

void VectorExtruder::first_fill(Polymer &poly) {
  bool is_overl = false;

  std::vector<Extruder> tmp_extruder;
  for (int i = 0; i < m_n_max_extr; ++i) {
    m_extr.place_extruder(poly);
    if( tmp_extruder.size() != 0){
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
  
  else{
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

void VectorExtruder::first_fill_field(Polymer &poly, FieldAction cohes_field_int) {
  bool is_overl = false;
  std::vector<Extruder> tmp_extruder;
  
  std::vector<int> subchain_cell;
  cohes_field_int.interaction();
  	
  for(auto &a : cohes_field_int.get_contact_cell()){
  // Quanti estrusori si possono attaccare in un solo step in una singola cella?
      int extr_per_step=1;
      for( int i = 0; i< extr_per_step; ++i){	  
        m_extr.place_extruder_cell(poly, cohes_field_int.monomer_min(a), cohes_field_int.monomer_max(a));
  	for (auto &j : tmp_extruder)
        	if (j.extr_overlap(m_extr)) {
         	 	is_overl = true;
          		break;
        	}
  	if (( extr_per_step*cohes_field_int.get_c(a.i, a.j, a.k)
	     *m_kon*integrator_timestep/poly.get_poly_nmonomers() ) > dist(mt) and
          !(is_overl))
        	tmp_extruder.push_back(m_extr);
    }
  }
  
  m_vector_extr.clear();
  for (const auto &i : tmp_extruder)
    m_vector_extr.push_back(std::make_unique<Extruder>(i));
}

void VectorExtruder::update_field(Polymer &poly, FieldAction cohes_field_int) {
  bool is_overl = false;
  std::vector<Extruder> tmp_extruder;
  for (const auto &i : m_vector_extr)
    // fill tmp_extruder only with extruder
    // that are not undbind
    if ((m_vector_extr.size() * m_koff * integrator_timestep) < dist(mt))
      tmp_extruder.push_back(*i);
  
  std::vector<int> subchain_cell;
  cohes_field_int.interaction();
  
   for(auto &a : cohes_field_int.get_contact_cell()){
  // Quanti estrusori si possono attaccare in un solo step in una singola cella?
      int extr_per_step=10;
      //sistemare: ipoteticamente posso attaccare infiniti estrusori in uno step
      for( int i = 0; i< extr_per_step - tmp_extruder.size(); ++i){	  
        m_extr.place_extruder_cell(poly, cohes_field_int.monomer_min(a), cohes_field_int.monomer_max(a));
  	for (auto &j : tmp_extruder)
        	if (j.extr_overlap(m_extr)) {
         	 	is_overl = true;
          		break;
        	}
  	if (( extr_per_step*cohes_field_int.get_c(a.i, a.j, a.k)
	     *m_kon*integrator_timestep/poly.get_poly_nmonomers() ) > dist(mt) and
          !(is_overl))
        	tmp_extruder.push_back(m_extr);
    }
  }
  
  m_vector_extr.clear();
  for (const auto &i : tmp_extruder)
    m_vector_extr.push_back(std::make_unique<Extruder>(i));
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
