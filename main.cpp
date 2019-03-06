/*
 * main.cpp
 *
 *  Created on: Feb 26, 2019
 *      Author: martanit
 */

#include "parameters.h"
#include "polymer.h"

int main() {

//	Parameters parm("data.dat");
	Parameters parm2;
	parm2.read_parm("data.dat");
	Polymer poly;	
	poly.poly_configuration();

	return 0;
}
