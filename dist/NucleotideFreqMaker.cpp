/*
 * BaseFreqMaker.cpp
 *
 *  Created on: Apr 24, 2012
 *      Author: Hani Zakaria Girgis, Ph.D. at the NCBI/NLM/NIH - USA
 */
#include "Util.h"

#include <string>
#include <vector>
#include <iostream>

using namespace std;

double count[] = { 0.0, 0.0, 0.0, 0.0 };

void readChromosome(string& chrom) {
	for (int i = 0; i < chrom.length(); i++) {
		switch (chrom[i]) {
		case 'A':
			count[0]++;
			break;
		case 'C':
			count[1]++;
			break;
		case 'G':
			count[2]++;
			break;
		case 'T':
			count[3]++;
			break;
		case 'a':
		case 'c':
		case 'g':
		case 't':
			cerr << "Error: nucleotides are expected to be large case." << endl;
			break;
		default:
			// Do nothing
			break;
		}
	}
}

void readGenome(string genomeDir) {
	vector<string> * chromList = Util::readChromList(genomeDir);

	vector<string> * infoList = new vector<string> ();
	infoList->reserve(chromList->size());

	vector<string> * seqList = new vector<string> ();
	seqList->reserve(chromList->size());

	for (int i = 0; i < chromList->size(); i++) {
		cerr << "Processing " << chromList->at(i) << " ..." << endl;
		Util::readFasta(chromList->at(i), infoList, seqList);
		Util::toUpperCase(seqList->at(i));
		readChromosome(seqList->at(i));
	}

	double sum = count[0] + count[1] + count[2] + count[3];
	cerr << "A: " << count[0] / sum << endl;
	cerr << "C: " << count[1] / sum << endl;
	cerr << "G: " << count[2] / sum << endl;
	cerr << "T: " << count[3] / sum << endl;
	cout << count[0] / sum << " " << count[1] / sum << " " << count[2] / sum
			<< " " << count[3] / sum << endl;

}

int main(int argc, char * argv[]) {
	if (argc != 2) {
		cerr << "Usage: " << argv[0] << " directory" << endl;
	} else {
		readGenome(argv[1]);
	}
	return 0;
}
