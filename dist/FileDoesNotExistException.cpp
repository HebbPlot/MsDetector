/*
 * FileDoesNotExistException.cpp
 *
 *  Created on: Apr 30, 2012
 *      Author: Hani Zakaria Girgis, Ph.D. at the NCBI/NLM/NIH - USA
 */

#include "FileDoesNotExistException.h"

#include <iostream>
#include <string>

using namespace std;

FileDoesNotExistException::FileDoesNotExistException(string massage) {
	cerr << massage << endl;
}

FileDoesNotExistException::~FileDoesNotExistException() {
	// TODO Auto-generated destructor stub
}
