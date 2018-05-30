/*
 * InvalidInputException.h
 *
 *  Created on: May 1, 2012
 *      Author: Hani Zakaria Girgis, Ph.D. at the NCBI/NLM/NIH - USA
 */

#ifndef INVALIDINPUTEXCEPTION_H_
#define INVALIDINPUTEXCEPTION_H_

#include<string>

using namespace std;

class InvalidInputException {
public:
	InvalidInputException(string);
	~InvalidInputException();
};

#endif /* INVALIDINPUTEXCEPTION_H_ */
