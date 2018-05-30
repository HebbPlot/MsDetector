//============================================================================
// Name        : MsDetector.cpp
// Author      : Hani Zakaria Girgis, Ph.D. at the NCBI/NLM/NIH - USA
//============================================================================

#include "Coordinator.h"
#include "Util.h"

#include <HMMlib/hmm_table.hpp>
#include <HMMlib/hmm_vector.hpp>
#include <HMMlib/hmm.hpp>

#include <iostream>
#include <string>
#include <fstream>
#include <vector>

using namespace std;

HMM<double> * buildHmm() {
	int K = 2; // number of states
	int M = 7; // alphabet size

	boost::shared_ptr<HMMVector<double> > pi_ptr(new HMMVector<double> (K));
	boost::shared_ptr<HMMMatrix<double> > T_ptr(new HMMMatrix<double> (K, K));
	boost::shared_ptr<HMMMatrix<double> > E_ptr(new HMMMatrix<double> (M, K));

	// initial probabilities
	HMMVector<double> &pi = *pi_ptr;
	double pSum = 39244 + 427;
	pi(0) = 39244 / pSum;
	pi(1) = 427 / pSum;
	// cout << pi(0) << " " << pi(1) << endl;

	// transitions
	HMMMatrix<double> &T = *T_ptr;
	double sumS0 = 19596025 + 4728;
	double sumS1 = 4667 + 190409;
	T(0, 0) = 19596025 / sumS0;
	T(0, 1) = 4728 / sumS0;
	T(1, 0) = 4667 / sumS1;
	T(1, 1) = 190409 / sumS1;
	// cout << T(0, 0) << " " << T(0, 1) << endl;
	// cout << T(1, 0) << " " << T(1, 1) << endl;

	HMMMatrix<double> &E = *E_ptr;
	double sumO0 = 0 + 0 + 8427 + 2506348 + 11124967 + 5265515 + 720665;
	E(0, 0) = 0 / sumO0;
	E(1, 0) = 0 / sumO0;
	E(2, 0) = 8427 / sumO0;
	E(3, 0) = 2506348 / sumO0;
	E(4, 0) = 11124967 / sumO0;
	E(5, 0) = 5265515 / sumO0;
	E(6, 0) = 720665 / sumO0;

	// Need to update the new numbers
	double sumO1 = 0 + 0 + 17 + 1952 + 15220 + 43762 + 148627;
	E(0, 1) = 0 / sumO1;
	E(1, 1) = 0 / sumO1;
	E(2, 1) = 17 / sumO1;
	E(3, 1) = 1952 / sumO1;
	E(4, 1) = 15220 / sumO1;
	E(5, 1) = 43762 / sumO1;
	E(6, 1) = 148627 / sumO1;
	// cout << E(0, 0) << " " << E(1, 0) << " " << E(2, 0) << " " << E(3, 0)
	//		<< " " << E(4, 0) << " " << E(5, 0) << " " << E(6, 0) << endl;
	// cout << E(0, 1) << " " << E(1, 1) << " " << E(2, 1) << " " << E(3, 1)
	//		<< " " << E(4, 1) << " " << E(5, 1) << " " << E(6, 1) << endl;
	// cout << "Constructing HMM" << endl;

	HMM<double> * myHMM = new HMM<double> (pi_ptr, T_ptr, E_ptr);
	return myHMM;
}

vector<double> * buildFilter() {
	double m1 = 25.537;
	double m2 = 5.7773;
	double s1 = 35.802;
	double s2 = 0.23168;
	double b = 6.4758;
	double w1 = 21.631;
	double w2 = 2.7629;

	vector<double> * paramList = new vector<double> ();
	paramList->push_back(m1);
	paramList->push_back(m2);
	paramList->push_back(s1);
	paramList->push_back(s2);
	paramList->push_back(b);
	paramList->push_back(w1);
	paramList->push_back(w2);

	return paramList;
}

void coordinate(string seqFile, string maskedFile, string repeatFile) {
	// Read sequence file
	Util::checkFile(seqFile);
	vector<string> * infoList = new vector<string> ();
	infoList->reserve(1000);
	vector<string> * seqList = new vector<string> ();
	seqList->reserve(1000);
	Util::readFasta(seqFile, infoList, seqList);

	// Construct HMM
	HMM<double> * hmmPtr = buildHmm();

	// Construct Filter
	vector<double> * filterPtr = buildFilter();

	// Parameters of MsDetector
	bool canShuffle = false;
	int motifLen = 6;
	int windowFactor = 4;
	string matrix = "Id";

	// Start with new output files
	Util::deleteFile(maskedFile);
	Util::deleteFile(repeatFile);

	// Process the input sequence(s)
	for (int i = 0; i < infoList->size(); i++) {
		Coordinator * coor = new Coordinator(seqList->at(i), infoList->at(i),
				motifLen, windowFactor, canShuffle);

		// Operations
		// 1: scoring
		coor->score(matrix);
		// 2: detection
		coor->detect(hmmPtr);
		// 3: filtering
		coor->filter(filterPtr);
		// 4-a: masking repeats
		coor->maskRepeats(maskedFile);
		// 4-b: writing repeats coordinates
		coor->writeRepeatsCoordinates(repeatFile);

		delete coor;
	}

	// Free memory
	delete hmmPtr;
	filterPtr->clear();
	delete filterPtr;
	infoList->clear();
	delete infoList;
	seqList->clear();
	delete seqList;
}

void printMassage(){
        cout << endl;
        cout << "MsDetector is designed and developed by Hani Zakaria Girgis under the supervision" << endl;
        cout << "of Sergey Sheetlin at the Spouge Research Group, the National Center for" << endl;
        cout << "Biotechnology Information, The National Institutes of Health, USA." << endl;
        cout << endl;
        cout << "The library HMMlib is utilized in MsDetector. HMMlib is covered by the GNU Lesser" << endl;
        cout << "General Public License." << endl;
        cout << endl;
        cout << "Version: 1.2" << endl;
        cout << endl;
}

int main(int argc, char * argv[]) {
	printMassage();
	// Usage massage
	string massage = string(
			"Usage: MsDetectorOptimized sequenceFile maskedFile msFile\n");
	massage.append("\tsequenceFile: input - sequence file\n");
	massage.append("\tmaskedFile: output - masked sequence file\n");
	massage.append("\tmsFile: output - MSs file\n");

	if (argc == 4) {
		coordinate(argv[1], argv[2], argv[3]);
	} else {
		cerr << massage << endl;
	}

	return 0;
}
