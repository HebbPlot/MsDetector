/*
 * InvalidOrderOfOperationsException.cpp
 *
 *  Created on: Apr 26, 2012
 *      Author: Hani Zakaria Girgis, Ph.D. at the NCBI/NLM/NIH - USA
 */

#include "InvalidOrderOfOperationsException.h"

#include <string>
#include <iostream>

using namespace std;

InvalidOrderOfOperationsException::InvalidOrderOfOperationsException(string massage) {
	cerr << massage << endl;
}

InvalidOrderOfOperationsException::~InvalidOrderOfOperationsException() {
	// TODO Auto-generated destructor stub
}
