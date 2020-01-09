#include "extruder.h"

// Try to place extruder randomly on the polymer
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

      else if (tmp_extr_pos != poly.get_poly_nmonomers() and
               m_ctcf[tmp_extr_pos + 1] == 0 and
               m_coupling_prob[tmp_extr_pos + 1] > tmp_coupling_try) {
        m_extruder_l = tmp_extr_pos;
        m_extruder_r = tmp_extr_pos + 1;
        set = 1;
      }
    }
  }
}

// Try to place extruder randomly on the polymer segment
// that is in a specific cell. The only information I need is the index
// of the monomers in the cell
void Extruder::place_extruder_cell(Polymer poly, int monomer_min, int monomer_max) {
  set = false;
  int p=0;
  while (!set) {
   std::cout << ++p << std::endl;
    try_extr_pos_cell = std::uniform_int_distribution<>(monomer_min, monomer_max);
    tmp_extr_pos = try_extr_pos_cell(mt);
    tmp_coupling_try = coupling_try(mt);

    if ((m_ctcf[tmp_extr_pos] == 0 and
         (m_ctcf[tmp_extr_pos - 1] == 0 or m_ctcf[tmp_extr_pos + 1] == 0)) and
        m_coupling_prob[tmp_extr_pos] > tmp_coupling_try) {

      if (tmp_extr_pos != monomer_min and m_ctcf[tmp_extr_pos - 1] == 0 and
          m_coupling_prob[tmp_extr_pos - 1] > tmp_coupling_try) {
        m_extruder_r = tmp_extr_pos;
        m_extruder_l = tmp_extr_pos - 1;
        set = true;
      }

      else if (tmp_extr_pos != monomer_max and
               m_ctcf[tmp_extr_pos + 1] == 0 and
               m_coupling_prob[tmp_extr_pos + 1] > tmp_coupling_try) {
        m_extruder_l = tmp_extr_pos;
        m_extruder_r = tmp_extr_pos + 1;
        set = true;
      }
    }
  }
}
bool Extruder::can_place_extr(Polymer poly, int monomer_min, int monomer_max){
  if(monomer_min == monomer_max)
	  return false;
  else {
	bool free=false;
  	for(int i=(monomer_min+1); i<monomer_max; ++i)
		  if(m_ctcf[i] == 0 and (m_ctcf[i+1]==0 or m_ctcf[i-1]==0)) free=true;
  	
	if(free==false) return false;
	else return true;
	}
}


bool Extruder::extr_overlap(Extruder &extr) {
  // two extruder cannot bind to monomer
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
