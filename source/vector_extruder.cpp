#include "vector_extruder.h"

void VectorExtruder::first_fill(Polymer &poly) {
  bool is_overl = false;

  std::vector<Extruder> tmp_extruder;
  for (int i = 0; i < m_n_max_extr; ++i) {
    m_extr.place_extruder(poly);
    // Fill m_vector extr with almost one extruder
    if (i == 0)
      tmp_extruder.push_back(m_extr);
    else {
      for (auto &j : tmp_extruder)
        if (j.extr_overlap(m_extr)) {
          is_overl = true;
          break;
        }
      if ((m_n_max_extr * m_kon * integrator_timestep /
           poly.get_poly_sphere()) > dist(mt) and
          !(is_overl))
        tmp_extruder.push_back(m_extr);
    }
  }
  m_vector_extr.clear();
  for (const auto &i : tmp_extruder)
    m_vector_extr.push_back(std::make_unique<Extruder>(i));
}

void VectorExtruder::update(Polymer &poly) {
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
    if ((m_n_max_extr * m_kon * integrator_timestep / poly.get_poly_sphere()) >
            dist(mt) and
        !(is_overl))
      tmp_extruder.push_back(m_extr);
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
