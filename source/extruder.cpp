#include "extruder.h"

void Extruder::place_extruder(Polymer poly){
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<> try_extruder_pos(1, poly.get_poly_sphere()-2);
    bool set = 0;
    int tmp_extr_pos;
    
    while(!set){
        tmp_extr_pos = try_extruder_pos(mt);
        if((m_ctcf[tmp_extr_pos]==0 and 
          (m_ctcf[tmp_extr_pos-1]==0 or  m_ctcf[tmp_extr_pos+1]==0)) 
            and m_coupling_prob[tmp_extr_pos]>0.5){
          
              if ( tmp_extr_pos!=0 and
                   m_ctcf[tmp_extr_pos-1]==0 and 
                   m_coupling_prob[tmp_extr_pos-1]>0.5 ) 
                {
                m_extruder_r = tmp_extr_pos;
                m_extruder_l = tmp_extr_pos-1;
                set = 1;
              }

              else if (tmp_extr_pos!=poly.get_poly_sphere() and 
                       m_ctcf[tmp_extr_pos+1]==0 and 
                       m_coupling_prob[tmp_extr_pos+1]>0.5 ){
                m_extruder_l = tmp_extr_pos;
                m_extruder_r = tmp_extr_pos+1;
                set = 1;
              }
        }
    }
}

// Input/Output function for write extruder positiom
bool print_r(Polymer & poly, Extruder& extr, std::string out_r )
{	
  // open stream to write xyz file
  std::ofstream output;
  output.open( out_r, std::fstream::app );

	// return error if read file fail
	if( output.fail() ) {
		throw "ERROR: Impossible to write xyz file "+out_r;
		return 1;
	}

/*		output << "\tAu\t\t\t" << poly.get_x(extr.m_extruder_l) 
                   << "\t\t\t" << poly.get_y(extr.m_extruder_l) 
                   << "\t\t\t" << poly.get_z(extr.m_extruder_l) 
                   << "\t\t\t" << poly.get_x(extr.m_extruder_r) 
                   << "\t\t\t" << poly.get_y(extr.m_extruder_r) 
                   << "\t\t\t" << poly.get_z(extr.m_extruder_r) << std::endl;
*/
    output << extr.m_extruder_l << "   " <<extr.m_extruder_r << std::endl;
  output.close();
  return 0;
}
