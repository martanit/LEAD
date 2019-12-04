/*
 * cohesin_field.h
 *
 *  Created on: December 1, 2019
 *  	Author: martanit
 */

#ifndef COHESIN_FIELD_H_
#define COHESIN_FIELD_H_

#include "parameters.h"

#include <cmath>
#include <fstream>
#include <random>
#include <vector>
#include <algorithm>
#include <memory>
class CohesinField {

public:
  // default constructor
  CohesinField(){};
  // assign passed parameters and construct the field
  CohesinField(Parameters);

  ~CohesinField();
  
  void operator=(const CohesinField &lhs) {
   	  
    for (auto &i : (lhs.m_cohesin_c)){
	m_cohesin_c[i].clear()    
      	m_cohesin_c[i][j][k] = std::make_unique<CohesinField>(*k);
  }
  
  // place randomly first cohesin configuration in filed
  void set_cohesin();
  const double &get_c() const { return c; }
  const double &get_field_step() const { return m_field_step; }
  const double &get_k_diff() const { return m_k_diff; }

private:
  
  // cohesin field parameters
  int m_poly_sphere = 100;   // number of interacting sphere
  double m_field_step;
  double m_k_diff;
  // cohesin that diffuse
  double c = 0.1;
   
  std::vector<std::vector<std::vector<double>>> m_cohesin_c;
  std::vector<std::vector<std::vector<double>>> m_cohesin_c_new;

  std::mt19937 mt{std::random_device{}()};
  std::uniform_real_distribution<double> uniform01;
  std::uniform_int_distribution<int> uniform05;

};

#endif /* COHESIN_FIELD_H_ */
