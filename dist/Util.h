/*
 * Util.h
 *
 *  Created on: Apr 24, 2012
 *      Author: Hani Zakaria Girgis, Ph.D. at NCBI/NLM/NIH - USA
 */

#ifndef UTIL_H_
#define UTIL_H_

#include <vector>
#include <string>


using namespace std;

class Util {
private:
	Util();
	~Util();

public:
	static string fileSeparator;
	static void readFasta(string, vector<string> *, vector<string> *);
	static vector<string> * readChromList(string);
	static void toUpperCase(string*);
	static void toUpperCase(string&);
	static string int2string(int);
	static string double2string(double);
	static void deleteFile(string);
	static void checkFile(string);
};

#endif /* UTIL_H_ */
