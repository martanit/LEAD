#ifndef UTILS_H_
#define UTILS_H_

#include <cmath>
// algorithm for periodic boundary conditions with side L=box
//static const double pbc(double r) { return  r - this->m_box * std::rint(r/this->m_box); };
static const double pbc(double r) { return  r; };

#endif /* UTILS_H_ */
