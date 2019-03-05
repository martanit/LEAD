/*
 * main.cpp
 *
 *  Created on: Feb 26, 2019
 *      Author: martanit
 */

#include "parameters.h"

int main() {

//	Parameters parm("data.dat");
	Parameters parm2;
	parm2.read_parm("data.dat");

	std::cout << parm2.get_temp() << std::endl;

	return 0;
}
