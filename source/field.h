/*
 * field.h
 *
 *  Created on: December 1, 2019
 *    Author: martanit
 *
 * field_action.h
 *
 *  Created on: Gen 03, 2020
 *    Author: martanit
 */

#ifndef FIELD_H_
#define FIELD_H_

#include "parameters.h"
#include "polymer.h"

#include <fstream>
#include <random>
#include <vector>

struct Cell;
class Field;
class FieldAction;

struct Cell
{
    int i, j, k;
};

class Field {

public:

    // default constructor
    Field() {};
    
    // assign passed parameters and construct the field
    Field(Parameters);
    
    ~Field() {};

    void operator=(const Field &rhs)
    {
        m_extruder_c = rhs.m_extruder_c;
        m_field_length = rhs.m_field_length;
        m_field_step = rhs.m_field_step;
        m_delta_c = rhs.m_delta_c;
        m_k_diff = rhs.m_k_diff;
    }

    // place randomly first extruders configuration in filed
    void first_extruders_field();
    
    // add a packet of extruders, aka a quantity of cohesin
    void add_delta_c(int, int, int);
    void sub_delta_c(int, int, int);

    //access funciton
    const double &get_init_c() const {
        return m_init_c;
    }
    const double &get_delta_c() const {
        return m_delta_c;
    }
    const double &get_c(int i, int j, int k) const {
        return m_extruder_c[i][j][k];
    }
    const double &get_field_step() const {
        return m_field_step;
    }
    const double &get_field_length() const {
        return m_field_length;
    }
    const double &get_k_diff() const {
        return m_k_diff;
    }

private:

    // extruder field parameters
    double m_field_length = 10.;
    double m_field_step = 1.;
    // initial extruders per cell
    double m_init_c = 46.;
    // extruders quantity that diffuse
    double m_delta_c = 1.;
    // extruder diffusion rate
    double m_k_diff = 1E-10;

    // Lattice of extruders
    std::vector<std::vector<std::vector<double>>> m_extruder_c;

    std::mt19937 mt{std::random_device{}()};
    std::uniform_real_distribution<double> uniform01;
    std::uniform_int_distribution<int> uniform05;

};

class FieldAction : public Field {

public:

    FieldAction() {};
    
    FieldAction(Parameters parm, Polymer poly) : Field(parm),
        m_poly(poly),
        scale_length(get_field_step()),
        box_length(get_field_length())
    {
        m_poly.set_cm();
        shift_x=m_poly.get_xcm();
        shift_y=m_poly.get_ycm();
        shift_z=m_poly.get_zcm();
        this->interaction();
    };
    
    ~FieldAction() {};

    void operator=(const FieldAction &rhs)
    {
        Field::operator=(rhs);
        m_contact_cell = rhs.m_contact_cell;
        m_poly = rhs.m_poly;
    }

    bool poly_in_cell(Cell);
    void interaction();
    void subchain_in_cell(Cell);

    int monomer_min(Cell);
    int monomer_max(Cell);

    // Access function
    const std::vector<Cell> &get_contact_cell() const {
        return m_contact_cell;
    };

    // Conversion field to space coordinates
    double shift_x, shift_y, shift_z;
    double scale_length;
    double box_length;

    const double x(const int i) const {
        return scale_length*i-get_field_length()/2.*scale_length+shift_x;
    };
    const double y(const int i) const {
        return scale_length*i-get_field_length()/2.*scale_length+shift_y;
    };
    const double z(const int i) const {
        return scale_length*i-get_field_length()/2.*scale_length+shift_z;
    };

    friend void print_field(FieldAction, std::string);

private:
    
    Polymer m_poly;
    std::vector<Cell> m_contact_cell;
    std::vector<int> m_poly_subchain;
};

#endif /*FIELD_H_*/
