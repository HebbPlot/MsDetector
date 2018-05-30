/*
 * FileDoesNotExistException.h
 *
 *  Created on: Apr 30, 2012
 *      Author: Hani Zakaria Girgis, Ph.D. at the NCBI/NLM/NIH - USA
 */

#ifndef FILEDOESNOTEXISTEXCEPTION_H_
#define FILEDOESNOTEXISTEXCEPTION_H_

#include <string>

using namespace std;

class FileDoesNotExistException {
public:
	FileDoesNotExistException(string);
	~FileDoesNotExistException();
};

#endif /* FILEDOESNOTEXISTEXCEPTION_H_ */
