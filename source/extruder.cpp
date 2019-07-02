#include "extruder.h"

void Extruder::place_extruder(Polymer poly) {
  set = 0.;
  while (!set) {
    tmp_extr_pos = try_extruder_pos(mt);
    tmp_coupling_try = coupling_try(mt);
    if ((m_ctcf[tmp_extr_pos] == 0 and
         (m_ctcf[tmp_extr_pos - 1] == 0 or m_ctcf[tmp_extr_pos + 1] == 0)) and
        m_coupling_prob[tmp_extr_pos] > tmp_coupling_try) {

      if (tmp_extr_pos != 0 and m_ctcf[tmp_extr_pos - 1] == 0 and
          m_coupling_prob[tmp_extr_pos - 1] > tmp_coupling_try) {
        m_extruder_r = tmp_extr_pos;
        m_extruder_l = tmp_extr_pos - 1;
        set = 1;
      }

      else if (tmp_extr_pos != poly.get_poly_sphere() and
               m_ctcf[tmp_extr_pos + 1] == 0 and
               m_coupling_prob[tmp_extr_pos + 1] > tmp_coupling_try) {
        m_extruder_l = tmp_extr_pos;
        m_extruder_r = tmp_extr_pos + 1;
        set = 1;
      }
    }
  }
}

bool Extruder::extr_overlap(Extruder &extr) {
  // two extruder cannot bind to sphere
  // that is already taken from another extr
  if ((this->m_extruder_l == extr.m_extruder_l) and
      (this->m_extruder_r == extr.m_extruder_r))
    return true;

  else
    return false;
}

// Input/Output function for write extruder positiom
bool print_r(Polymer &poly, Extruder &extr, std::string out_r) {
  // open stream to write xyz file
  std::ofstream output;
  output.open(out_r, std::fstream::app);

  // return error if read file fail
  if (output.fail()) {
    throw "ERROR: Impossible to write xyz file " + out_r;
    return 1;
  }
  output << extr.m_extruder_l << "   " << extr.m_extruder_r << std::endl;
  output.close();
  return 0;
}
