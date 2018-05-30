/*
 * Coordinator.h
 *
 *  Created on: Oct 18, 2011
 *      Author: Hani Zakaria Girgis, Ph.D. at the NCBI/NLM/NIH - USA
 */
#ifndef COORDINATOR_H_
#define COORDINATOR_H_

#include <string>
#include <vector>
#include <HMMlib/hmm.hpp>

#include "ScoresMaker.h"
#include "CompositionMatrix.h"

using namespace std;
using namespace hmmlib;

class Coordinator {
public:
	Coordinator(string, string, int, int, bool);
	virtual ~Coordinator();
	void score(string);
	void score(CompositionMatrix *);
	void writeScores(string);
	void detect(HMM<double> *);
	void filter(vector<double> *);
	void filter(vector<double> *, double);
	void maskRepeats(string);
	void writeRepeatsCoordinates(string);

private:
	// string seqFile;
	int length;
	int winLength;
	string info;
	string base;
	const char * sequence; // Start of the base
	int gap;
	string chrm;
	int coorStart;
	int coorEnd;
	int step;
	int minSegLength;
	bool canDetect;
	bool canFilter;
	bool canOutput;
	vector<vector<int> *> * segment;
	vector<ScoresMaker *> * scorers;
	vector<vector<unsigned int> *> * repeats;
	vector<vector<double> *> * repIndex;
	HMM<double> * myHMM;

	void removeN();
	void makeSegments();
	void shuffle();
	void detectRepeats(ScoresMaker *);
	void extractRepeatIndex(ScoresMaker *, int);
	double getMean(int *, int, int);
	double getSum(int *, int, int);
};

#endif /* COORDINATOR_H_ */
