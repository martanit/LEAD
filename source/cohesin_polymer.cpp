#include "cohesin_polymer.h"

void CohesinPolymer::poly_field_interaction() {
		bool is_poly_in_cell = false;
		Cell a;
		m_contact_cell.clear();
		for(int i = 0; i< box_length; ++i)
				for(int j = 0; j< box_length; ++j)
					for(int k = 0; k< box_length; ++k){
						a = {i,j,k};
						is_poly_in_cell = poly_in_cell({i,j,k});
					}		
						if(is_poly_in_cell){
								m_contact_cell.push_back(a);
						}	
}

//void CohesinPolymer::update_poly_field_int(){
// Assumo che il poly si muova poco e
// controllo solo i primi vicini, più efficente :D!
// da IMPLEMENTARE!!

//}

bool CohesinPolymer::poly_in_cell(Cell a){
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

std::vector<int> CohesinPolymer::subchain_in_cell(Cell a){
		std::vector<int> poly_subchain;
	  for(int l = 0; l<m_poly.get_poly_nmonomers(); ++l)
				if((x(a.i) < m_poly.get_x(l) and m_poly.get_x(l) < x(a.i)+scale_length) and
				  (y(a.j) < m_poly.get_y(l) and m_poly.get_y(l) < y(a.j)+scale_length) and
  		    (z(a.k) < m_poly.get_z(l) and m_poly.get_z(l) < z(a.k)+scale_length))
						poly_subchain.push_back(l);
  return poly_subchain;
}

const int &CohesinPolymer::monomer_min(const Cell a)  {
	return *std::min_element(this->subchain_in_cell(a).begin(), this->subchain_in_cell(a).end());
} 

const int &CohesinPolymer::monomer_max(const Cell a)  {
	return *std::max_element(this->subchain_in_cell(a).begin(), this->subchain_in_cell(a).end());
}

