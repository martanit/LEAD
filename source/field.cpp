#include "field.h"

Field::Field(Parameters parm)
    :m_field_step(parm.get_field_step()),
     m_field_length(parm.get_field_length()),
     m_Dextr_free(parm.get_Dextr_free()),
     uniform01(0., 1.),
     uniform05(0, 5),
     m_extruder_c(m_field_length, std::vector<std::vector<double>>(m_field_length, std::vector<double>(m_field_length))) {
    this->first_extruders_field();
}

// assuming uniform distribution for extruder concentration
void Field::first_extruders_field() {
    for(auto &i : m_extruder_c)
        for(auto &j : i)
            for(auto &k : j)
                k = m_init_c;
}

void Field::add_delta_c(int i, int j, int k) {
    m_extruder_c[i][j][k] += m_delta_c;
}

void Field::sub_delta_c(int i, int j, int k) {
    m_extruder_c[i][j][k] -= m_delta_c;
}

bool FieldAction::poly_in_cell(Cell a) {
    bool in_cell=false;
    for(int l=0; l<m_poly.get_poly_nmonomers(); ++l) {
        if((x(a.i) < m_poly.get_x(l) and
                m_poly.get_x(l) < (x(a.i)+scale_length)) and
                (y(a.j) < m_poly.get_y(l) and
                 m_poly.get_y(l) < (y(a.j)+scale_length)) and
                (z(a.k) < m_poly.get_z(l) and
                 m_poly.get_z(l) < (z(a.k)+scale_length))) {
            in_cell = true;
            break;
        }
    }
    return in_cell;
}

void FieldAction::interaction() {
    Cell a;
    m_contact_cell.clear();
    for(int i = 0; i< box_length; ++i)
        for(int j = 0; j< box_length; ++j)
            for(int k = 0; k< box_length; ++k) {
                a = {i,j,k};
                if(poly_in_cell(a))
                    m_contact_cell.push_back(a);
            }
}

//void FieldAction::update_poly_field_int(){
// Assumo che il poly si muova poco e
// controllo solo i primi vicini, più efficente :D!
// da IMPLEMENTARE!!

//}

void FieldAction::subchain_in_cell(Cell a) {
    m_poly_subchain.clear();
    for(int l = 0; l<m_poly.get_poly_nmonomers(); ++l)
        if((x(a.i) < m_poly.get_x(l) and m_poly.get_x(l) < x(a.i)+scale_length) and
                (y(a.j) < m_poly.get_y(l) and m_poly.get_y(l) < y(a.j)+scale_length) and
                (z(a.k) < m_poly.get_z(l) and m_poly.get_z(l) < z(a.k)+scale_length))
            m_poly_subchain.push_back(l);
}


int FieldAction::monomer_min(Cell a)  {
    this->subchain_in_cell(a);
    if(m_poly_subchain.size() == 1) return m_poly_subchain.at(0);
    else return *std::min_element(std::begin(m_poly_subchain), std::end(m_poly_subchain));
}

int FieldAction::monomer_max(Cell a)  {
    this->subchain_in_cell(a);
    if(m_poly_subchain.size() == 1) return m_poly_subchain.at(0);
    else return *std::max_element(std::begin(m_poly_subchain), std::end(m_poly_subchain));
}

void print_field(FieldAction field, std::string out_field) {
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
                if(field.get_c(i,j,k)>=46)
                    output << "\tNan" << "\t\t\t"
                           << field.x(i) << "\t\t\t"
                           << field.y(j) << "\t\t\t"
                           << field.z(k) << std::endl;
                else
                    output << "\tAu" << "\t\t\t"
                           << field.x(i) << "\t\t\t"
                           << field.y(j) << "\t\t\t"
                           << field.z(k) << std::endl;
    output.close();
}
