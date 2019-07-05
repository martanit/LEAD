/*
 * potential.h
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

class Polymer {

public:
  // assign passed parameters to polymer
  Polymer(Parameters);
  // contruct polymer from x, y, z coordinate file and assign parameters
  Polymer(Parameters, std::string);

  ~Polymer();

  void set_size();
  // place first sphere of polymer
  void first_sphere();
  // place randomly other sphere
  void poly_configuration();
  // check if spheres overlap (self avoiding polymer)
  bool is_overlap(int);

  // calculate distance between atoms
  double dist(const int &, const int &);

  // fuction to access sphere coordinates, velocities and forces
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

  void reset_force();
  void add_force(const int &, const double &, const double &, const double &);

  void reset_energy();
  void add_energy(const double &);

  // function to access polymer parameters
  const int &get_poly_sphere() const { return m_poly_sphere; };
  const float &get_poly_dist() const { return m_poly_dist; };
  const float &get_spring() const { return m_poly_spring; };

  friend bool print_xyz(Polymer &, std::string);
  friend bool read_xyz(Polymer &, std::string);

private:
  double d = 0.;

  // polymer parameters
  int m_poly_sphere = 100;   // number of interacting sphere
  float m_poly_dist = 1.;    // length of spring [50*nm]
  float m_poly_spring = 100.; // spring constant of harmonic oscillator (spring)
  float m_poly_rinit = 0.5;  // hard core radius of sphere
  std::string m_poly_init = "init.dat";

  std::vector<double> m_poly_x, m_poly_y, m_poly_z;
  std::vector<double> m_poly_fx, m_poly_fy, m_poly_fz;

  double m_poly_e = 0.;

  std::mt19937 mt{std::random_device{}()};
  std::uniform_real_distribution<double> uniform01;
};

#endif /* POLYMER_H_ */
