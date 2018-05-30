/*
 * Util.cpp
 *
 *  Created on: Apr 24, 2012
 *      Author: Hani Zakaria Girgis, Ph.D. at NCBI/NLM/NIH - USA
 *      This class has a collection of utilities.
 */
#include "Util.h"
#include "FileDoesNotExistException.h"
#include "InvalidInputException.h"

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <dirent.h>

using namespace std;

Util::Util() {
	// TODO Auto-generated constructor stub

}

Util::~Util() {
	// TODO Auto-generated destructor stub
}

string Util::fileSeparator("/");

void Util::readFasta(string seqFile, vector<string> * infoList,
		vector<string> * seqList) {
	ifstream in(seqFile.c_str());
	string info;

	bool isFirst = true;
	string * basePtr = new string("");
	while (in.good()) {
		string line;
		getline(in, line);
		if (line[0] == '>') {
			infoList->push_back(line);
			if (!isFirst) {
				seqList->push_back(*basePtr);
				basePtr = new string("");
			} else {
				isFirst = false;
			}
		} else {
			basePtr->append(line);
		}
	}
	seqList->push_back(*basePtr);
	in.close();

	// Post condition
	if (infoList->size() != seqList->size()) {
		cerr << "Error while reading the fasta input file. "
				<< "Header count = " << infoList->size() << " "
				<< "Sequence count = " << seqList->size() << endl;
		exit(1);
	}
}

vector<string> * Util::readChromList(string genomeDir) {
	vector<string> * chromList = new vector<string> ();

	// This function may not be platform-independent
	// Credit: http://www.cplusplus.com/forum/beginner/9173/
	DIR * dirPtr = opendir(genomeDir.c_str());

	struct dirent * entry;
	while (entry = readdir(dirPtr)) {
		string file(entry->d_name);
		// Credit: http://stackoverflow.com/questions/51949/how-to-get-file-extension-from-string-in-c
		if (file.substr(file.find_last_of(".") + 1) == "fa") {
			chromList->push_back(genomeDir + fileSeparator + entry->d_name);
		}
	}

	closedir(dirPtr);

	return chromList;
}

// This method will modify the contents of its parameter basePtr!
void Util::toUpperCase(string * basePtr) {
	string base = *basePtr;
	// Convert alphabet to upper case
	for (int i = 0; i < base.length(); i++) {
		base[i] = toupper(base[i]);
	}
}

void Util::toUpperCase(string& base) {
	// Convert alphabet to upper case
	for (int i = 0; i < base.length(); i++) {
		base[i] = toupper(base[i]);
	}
}

// credit: http://stackoverflow.com/questions/228005/alternative-to-itoa-for-converting-integer-to-string-c
string Util::int2string(int i) {
	string s;
	stringstream out;
	out << i;
	s = out.str();
	return s;
}

// Need to use templates
string Util::double2string(double i) {
	string s;
	stringstream out;
	out << i;
	s = out.str();
	return s;
}

void Util::checkFile(string fileName) {
	ifstream f1(fileName.c_str());
	if (!f1) {
		string massage = string("ERROR: ");
		massage.append(fileName);
		massage.append(" does not exist.\n");
		throw FileDoesNotExistException(massage);
	}
	f1.close();
}

void Util::deleteFile(string fileName) {
	ifstream f1(fileName.c_str());
	if (f1) {
		if (remove(fileName.c_str()) != 0) {
			cerr << "Could not remove: " << fileName << endl;
		} else {
			// cout << "Deleting: " << fileName << endl;
		}
	}
	f1.close();
}
