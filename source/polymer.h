/*
 * polymer.h
 *
 *  Created on: March 6, 2019
 *  	Author: martanit
 */

#ifndef POLYMER_H_
#define POLYMER_H_

#include "parameters.h"

#include <cmath>
#include <fstream>
#include <random>
#include <vector>
#include <algorithm>

class Polymer {

public:
  // default constructor (called only for m_poly_old in dynamics.h)
  Polymer(){};
  // assign passed parameters to polymer
  Polymer(Parameters);
  // contruct polymer from x, y, z coordinate file and assign parameters
  Polymer(Parameters, std::string);

  ~Polymer();

  void set_size();
  // place first monomer of polymer
  void first_monomer();
  // place randomly other monomers
  void poly_configuration();
  // check if monomers overlap (self avoiding polymer)
  bool is_overlap(int);

  // calculate distance between monomers
  double dist(const int &, const int &);

  // fuction to access monomers coordinates, velocities and forces
  const double &get_x(int i) const { return m_poly_x[i]; };
  const double &get_y(int i) const { return m_poly_y[i]; };
  const double &get_z(int i) const { return m_poly_z[i]; };

  const double &get_fx(int i) const { return m_poly_fx[i]; };
  const double &get_fy(int i) const { return m_poly_fy[i]; };
  const double &get_fz(int i) const { return m_poly_fz[i]; };

  const double &get_energy() const { return m_poly_e; };

  void set_x(double x, int i) { m_poly_x[i] = x; };
  void set_y(double y, int i) { m_poly_y[i] = y; };
  void set_z(double z, int i) { m_poly_z[i] = z; };

  void center();
  void set_cm();

  const double &get_xcm() const { return x_cm; };
  const double &get_ycm() const { return y_cm; };
  const double &get_zcm() const { return z_cm; };

  void reset_force();
  void add_force(const int &, const double &, const double &, const double &);

  void reset_energy();
  void add_energy(const double &);

  // function to access polymer parameters
  const int &get_poly_nmonomers() const { return m_poly_nmonomers; };
  const float &get_poly_d() const { return m_poly_d; };
  const float &get_spring() const { return m_poly_spring; };

  friend bool print_xyz(Polymer &, std::string);
  friend bool read_xyz(Polymer &, std::string);

private:
  
  double d = 0.;
  double sum_x, sum_y, sum_z;
  double x_cm, y_cm, z_cm;

  // polymer parameters
  int m_poly_nmonomers = 100;   // number of interacting monomers
  float m_poly_d = 0.89;    // diameter of monomers [50*nm]
  float m_poly_spring = 100.; // spring constant of harmonic oscillator (spring)
  std::string m_poly_init = "init.dat";

  std::vector<double> m_poly_x, m_poly_y, m_poly_z;
  std::vector<double> m_poly_fx, m_poly_fy, m_poly_fz;

  double m_poly_e = 0.;

  std::mt19937 mt{std::random_device{}()};
  std::uniform_real_distribution<double> uniform01;
};

#endif /* POLYMER_H_ */
