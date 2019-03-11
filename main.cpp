/*
 * main.cpp
 *
 *  Created on: Feb 26, 2019
 *      Author: martanit
 */

#include "parameters.h"
#include "polymer.h"

int main() {
	
	Polymer poly2("data.dat");
	poly2.poly_configuration();
	poly2.print_xyz( "prova.xyz");
	return 0;
}
