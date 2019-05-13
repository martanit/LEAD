#include <cmath>
#include "parameters.h"

class Utils 
{  
  public: 
    Utils() {};
    Utils(Parameters parm) : m_box(parm.get_box()) { };  

    ~Utils(){ };
    // algorithm for periodic boundary conditions with side L=box
    // double pbc(double r) { return  r - m_box * std::rint(r/m_box); };
   double pbc(double r) { return  r; };
 
  private:
    double m_box = 10.;
};
