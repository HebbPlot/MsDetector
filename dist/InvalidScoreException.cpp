/*
 * InvalidScoreException.cpp
 *
 *  Created on: Apr 27, 2012
 *      Author: Hani Zakaria Girgis, Ph.D. at the NCBI/NLM/NIH - USA
 */

#include "InvalidScoreException.h"

#include <string>
#include <iostream>

using namespace std;

InvalidScoreException::InvalidScoreException(string massage) {
	cerr << massage << endl;
}

InvalidScoreException::~InvalidScoreException() {
	// TODO Auto-generated destructor stub
}
