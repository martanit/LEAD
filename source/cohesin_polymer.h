/*
 * cohesin_polymer.h
 *
 *  Created on: Gen 03, 2020
 *  	Author: martanit
 */
#ifndef COHESIN_POLYMER_H_
#define COHESIN_POLYMER_H_

#include <fstream>
#include <random>
#include <vector>

#include "extruder_field.h"
#include "polymer.h"

struct Cell {
	int i, j, k;
};

class CohesinPolymer : public ExtruderField {
public:
  CohesinPolymer(){};
  CohesinPolymer(Parameters parm, Polymer poly) : ExtruderField(parm)
  {
			m_poly = poly;
			scale_length = get_field_step();
			box_length = get_field_length();
			shift_x=1; 
			shift_y=1;
		       	shift_z=1;
			this->poly_field_interaction();
  };
  ~CohesinPolymer(){};

  void operator=(const CohesinPolymer &rhs) {
	ExtruderField::operator=(rhs);
  	m_contact_cell = rhs.m_contact_cell;
	m_poly = rhs.m_poly;
  }
  void poly_field_interaction();
  bool poly_in_cell(Cell);
  std::vector<int> subchain_in_cell(Cell);
  
  const int &monomer_min(const Cell);
  const int &monomer_max(const Cell);
	
  // Access function
  const std::vector<Cell> &get_contact_cell() const { return m_contact_cell; };

  // Conversion field to space coordinates
  double scale_length, shift_x=1, shift_y=1, shift_z=1;
  double box_length;
  const double x(const int i) const { return i*scale_length + shift_x; };
  const double y(const int i) const { return i*scale_length + shift_y; };
  const double z(const int i) const { return i*scale_length + shift_z; };

private:
  Polymer m_poly;
  std::vector<Cell> m_contact_cell;
		
};

#endif /*COHESIN_POLYMER_H_*/
