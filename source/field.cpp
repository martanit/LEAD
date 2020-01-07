#include "field.h"

Field::Field(Parameters parm)
      :m_field_step(parm.get_field_step()),
	m_field_length(parm.get_field_length()),
      m_k_diff(parm.get_k_diff()),
      uniform01(0., 1.),
      uniform05(0, 5),
      m_extruder_c(m_field_length, std::vector<std::vector<double>>(m_field_length, std::vector<double>(m_field_length))){
	this->first_extruders_field();
}

// assuming uniform distribution for extruder concentration
void Field::first_extruders_field() {
    for(auto &i : m_extruder_c)
        for(auto &j : i)
            for(auto &k : j)
                k = uniform01(mt);
}

void Field::add_delta_c(int i, int j, int k) {
	m_extruder_c[i][j][k] += m_delta_c;
}

void Field::sub_delta_c(int i, int j, int k) {
	m_extruder_c[i][j][k] -= m_delta_c;
}

void FieldAction::interaction() {
		bool is_poly_in_cell = false;
		Cell a;
		m_contact_cell.clear();
		for(int i = 0; i< box_length; ++i)
				for(int j = 0; j< box_length; ++j)
					for(int k = 0; k< box_length; ++k){
					  a = {i,j,k};
					  is_poly_in_cell = poly_in_cell(a);
					  if(is_poly_in_cell) break;
					}		
					if(is_poly_in_cell)
						m_contact_cell.push_back(a);
					
				//std::cout << m_contact_cell.size() << std::endl;
}

//void FieldAction::update_poly_field_int(){
// Assumo che il poly si muova poco e
// controllo solo i primi vicini, più efficente :D!
// da IMPLEMENTARE!!

//}

bool FieldAction::poly_in_cell(Cell a){
		for(int l=0; l<m_poly.get_poly_nmonomers(); ++l)
				if((x(a.i) < m_poly.get_x(l) and 
						m_poly.get_x(l) < x(a.i)+scale_length) and
				   (y(a.j) < m_poly.get_y(l) and 
						m_poly.get_y(l) < y(a.j)+scale_length) and
				   (z(a.k) < m_poly.get_z(l) and 
						m_poly.get_z(l) < z(a.k)+scale_length)){
						return true;
				}
		else return false;
}

std::vector<int> FieldAction::subchain_in_cell(Cell a){
		std::vector<int> poly_subchain;
	  for(int l = 0; l<m_poly.get_poly_nmonomers(); ++l)
				if((x(a.i) < m_poly.get_x(l) and m_poly.get_x(l) < x(a.i)+scale_length) and
				  (y(a.j) < m_poly.get_y(l) and m_poly.get_y(l) < y(a.j)+scale_length) and
  		    (z(a.k) < m_poly.get_z(l) and m_poly.get_z(l) < z(a.k)+scale_length))
						poly_subchain.push_back(l);
  return poly_subchain;
}

const int &FieldAction::monomer_min(const Cell a)  {
	return *std::min_element(this->subchain_in_cell(a).begin(), this->subchain_in_cell(a).end());
} 

const int &FieldAction::monomer_max(const Cell a)  {
	return *std::max_element(this->subchain_in_cell(a).begin(), this->subchain_in_cell(a).end());
}

void print_field(FieldAction field, std::string out_field){
  	std::ofstream output;
  	output.open(out_field, std::fstream::app);

  	// return error if read file fail
  	if (output.fail()) {
    	throw "ERROR: Impossible to write field xyz file " + out_field;
  	}
	
	output << std::pow(field.get_field_length(),3) << std::endl<<std::endl;
	for(int i=0; i<field.get_field_length(); ++i)
    	    	for(int j=0; j<field.get_field_length(); ++j)
    	    		for(int k=0; k<field.get_field_length(); ++k)
			if(field.get_c(i,j,k)>=0.1)
		 	   output << "\tAu" << "\t\t\t"
		      			<< field.x(i) << "\t\t\t" 
		       			<< field.y(j) << "\t\t\t"
		       			<< field.z(k) << std::endl;
			else    		
		 	   output << "\tNan" << "\t\t\t"
		      			<< field.x(i) << "\t\t\t" 
		       			<< field.y(j) << "\t\t\t"
		       			<< field.z(k) << std::endl;
  	output.close();
    }

