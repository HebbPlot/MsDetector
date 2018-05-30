/*
 * InvalidOrderOfOperationsException.h
 *
 *  Created on: Apr 26, 2012
 *      Author: Hani Zakaria Girgis, Ph.D. at the NCBI/NLM/NIH - USA
 */

#ifndef INVALIDORDEROFOPERATIONSEXCEPTION_H_
#define INVALIDORDEROFOPERATIONSEXCEPTION_H_

#include <string>

using namespace std;

class InvalidOrderOfOperationsException {
public:
	InvalidOrderOfOperationsException(string);
	~InvalidOrderOfOperationsException();
};

#endif /* INVALIDORDEROFOPERATIONSEXCEPTION_H_ */
