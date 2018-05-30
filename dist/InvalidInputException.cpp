/*
 * InvalidInputException.cpp
 *
 *  Created on: May 1, 2012
 *      Author: Hani Zakaria Girgis, Ph.D. at the NCBI/NLM/NIH - USA
 */

#include "InvalidInputException.h"

#include <string>
#include <iostream>

using namespace std;

InvalidInputException::InvalidInputException(string msg) {
	cerr << msg << endl;
}

InvalidInputException::~InvalidInputException() {
	// TODO Auto-generated destructor stub
}
