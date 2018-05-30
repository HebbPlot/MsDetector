/*
 * Coordinator.cpp
 *
 *  Created on: Oct 18, 2011
 *      Author: Hani Zakaria Girgis, Ph.D. at the NCBI/NLM/NIH - USA
 */
#include "Coordinator.h"

#include "ScoresMaker.h"
#include "ScoresMatrixMaker.h"
#include "ScoresBiasMaker.h"
#include "CompositionMatrix.h"

#include "InvalidOrderOfOperationsException.h"

#include <HMMlib/hmm_table.hpp>
#include <HMMlib/hmm_vector.hpp>
#include <HMMlib/hmm.hpp>

#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <map>
#include <ctype.h>
#include <math.h>
#include <time.h>

using namespace std;
using namespace hmmlib;

/**
 * seq: nucleotide sequence
 * info: header of sequence
 * l: word length
 * factor: number of words in a half window
 * mtrx: the name of the matrix to be used to score the sequence
 */
Coordinator::Coordinator(string seq, string infoIn, int l, int factor,
		bool canShuffle) {
	base = seq;
	info = infoIn;
	length = l;
	// winLength = length * 4;
	winLength = length * factor;
	gap = 1;
	minSegLength = 20;
	step = 500;
	canDetect = false;
	canFilter = false;
	canOutput = false;

	// Break down the sequence header
	int colIndex = info.find_first_of(':');
	int dashIndex = info.find_first_of('-');
	if (colIndex < 0 || dashIndex < 0) {
		chrm = info.substr(1, info.length() - 1);
		coorStart = 0;
		coorEnd = seq.length();
	} else {
		chrm = info.substr(1, colIndex - 1);
		coorStart = atoi(
				info.substr(colIndex + 1, dashIndex - colIndex - 1).c_str());
		coorEnd = atoi(info.substr(dashIndex + 1).c_str());
	}

	segment = new vector<vector<int> *> ();
	segment->reserve(100);

	removeN();

	if (canShuffle) {
		shuffle();
	}

	makeSegments();

	repeats = new vector<vector<unsigned int> *> ();
	repeats->reserve(segment->size());

	repIndex = new vector<vector<double> *> ();
	repIndex->reserve(1000);

	sequence = base.c_str();

	scorers = new vector<ScoresMaker *> ();
	scorers->reserve(segment->size());
}

Coordinator::~Coordinator() {
	// cout << "Deleting memory used by the coordinator object .." << endl;

	segment->clear();
	delete segment;

	scorers->clear();
	delete scorers;

	repeats->clear();
	delete repeats;

	for (int i = 0; i < repIndex->size(); i++) {
		repIndex->at(i)->clear();
	}
	delete repIndex;
}

// Segment coordinates are inclusive [s,e]
// Segment consisting of less than 20 non-N bases are discarded
void Coordinator::removeN() {
	// Convert alphabet to upper case
	for (int i = 0; i < base.length(); i++) {
		base[i] = toupper(base[i]);
	}

	// Store non-N index
	int start = -1;
	for (int i = 0; i < base.size(); i++) {
		if (base[i] != 'N' && start == -1) {
			start = i;
		} else if (base[i] == 'N' && start != -1) {
			// Check the length of the segment before inserting
			if (i - start >= minSegLength) {
				vector<int> * v = new vector<int> ();
				v->push_back(start);
				v->push_back(i - 1);
				segment->push_back(v);
				start = -1;
			}
		} else if (i == base.size() - 1 && base[i] != 'N' && start != -1) {
			// Check the length of the segment before inserting
			if (i - start + 1 >= minSegLength) {
				vector<int> * v = new vector<int> ();
				v->push_back(start);
				v->push_back(i);
				segment->push_back(v);
				start = -1;
			}
		}
	}
}

// This method divides each segment into 1-mbp fragments.
// The last fragment can be longer than 1 mbp.
void Coordinator::makeSegments() {
	int step = 1000000;

	vector<vector<int> *> * temp = new vector<vector<int> *> ();
	for (int i = 0; i < segment->size(); i++) {
		int s = segment->at(i)->at(0);
		int e = segment->at(i)->at(1);
		if (e - s + 1 > step) {
			int fragNum = (int) (e - s + 1) / step;

			for (int h = 0; h < fragNum; h++) {
				int fragStart = s + (h * step);
				int fragEnd = (h == fragNum - 1) ? e : fragStart + step - 1;
				vector<int> * v = new vector<int> ();
				v->push_back(fragStart);
				v->push_back(fragEnd);
				temp->push_back(v);
			}
		} else if (e - s + 1 <= length) {
			// This a very short segment. Do nothing.
		} else {
			vector<int> * v = new vector<int> ();
			v->push_back(segment->at(i)->at(0));
			v->push_back(segment->at(i)->at(1));
			temp->push_back(v);
		}
	}

	segment->clear();
	segment = temp;
}

/**
 * Based on code found on the following link:
 * http://www.fredosaurus.com/notes-cpp/misc/random-shuffle.html
 **/
void Coordinator::shuffle() {
	for (int j = 0; j < segment->size(); j++) {
		int s = segment->at(j)->at(0);
		int e = segment->at(j)->at(1);
		int l = e - s + 1;

		srand(time(NULL));

		for (int i = 0; i < l - 1; i++) {
			int r = s + i + (rand() % (l - i)); // Random remaining position.

			int temp = base[s + i];
			base[s + i] = base[r];
			base[r] = temp;
		}
	}
}

/**
 * Handles the Id and the Trans matrixes
 */
void Coordinator::score(string matrix) {
	for (int i = 0; i < segment->size(); i++) {
		ScoresMaker * scoresMaker;
		if (matrix.compare("Id") == 0) {
			scoresMaker = new ScoresMaker(sequence, segment->at(i)->at(0),
					segment->at(i)->at(1), length, winLength);
		} else if (matrix.compare("Trans") == 0) {
			scoresMaker = new ScoresMatrixMaker(sequence,
					segment->at(i)->at(0), segment->at(i)->at(1), length,
					winLength);
		}

		scoresMaker->start();
		scorers->push_back(scoresMaker);
	}
	// Can detect now
	canDetect = true;
}

/**
 * Handles the Comp and the TransComp matrixes
 */
void Coordinator::score(CompositionMatrix * mtrx) {
	for (int i = 0; i < segment->size(); i++) {
		ScoresMaker * scoresMaker = new ScoresBiasMaker(sequence,
				segment->at(i)->at(0), segment->at(i)->at(1), length,
				winLength, mtrx);

		scoresMaker->start();
		scorers->push_back(scoresMaker);
	}
	// Can detect now
	canDetect = true;
}

void Coordinator::detect(HMM<double> * hmmPtr) {
	if (!canDetect) {
		throw InvalidOrderOfOperationsException(
				"Detection must occur after scoring.");
	}

	myHMM = hmmPtr;
	for (int i = 0; i < segment->size(); i++) {
		ScoresMaker * scoresMaker = scorers->at(i);
		detectRepeats(scoresMaker);
		extractRepeatIndex(scoresMaker, i);
	}

	// Can  now
	canFilter = true;
	canOutput = true;
}

void Coordinator::detectRepeats(ScoresMaker * scoresMaker) {
	int * scores = scoresMaker->getScores();
	int size = scoresMaker->getScoresSize();

	vector<unsigned int> * states = new vector<unsigned int> ();
	for (int s = 0; s < size; s += step) {
		int e = (s + step > size) ? size : s + step;

		// Index system [s, e[
		vector<unsigned int> * obs_ptr = new vector<unsigned int> (&scores[s],
				&scores[e]);
		vector<unsigned int> * hiddenseq_ptr = new vector<unsigned int> (e - s);

		myHMM->viterbi(*obs_ptr, *hiddenseq_ptr);

		states->insert(states->end(), hiddenseq_ptr->begin(),
				hiddenseq_ptr->end());

		obs_ptr->clear();
		delete obs_ptr;

		hiddenseq_ptr->clear();
		delete hiddenseq_ptr;
	}

	repeats->push_back(states);
}

void Coordinator::writeScores(string scoresFile) {
	ofstream out;
	out.open(scoresFile.c_str(), ios::out | ios::app);

	if (out) {
		for (int i = 0; i < scorers->size(); i++) {
			ScoresMaker * scoresMaker = scorers->at(i);

			int * scores = scoresMaker->getScores();
			int size = scoresMaker->getScoresSize();
			int start = scoresMaker->getStart();

			for (int s = 0; s < size; s += step) {
				int e = (s + step > size) ? size : s + step;
				out << chrm << ":" << start + s << "-" << start + e << " ";
				for (int j = s; j < e - 1; j++) {
					out << scores[j] << ",";
				}
				out << scores[e - 1] << endl;
			}
		}
		out.close();
	} else {
		cout << "ERROR: " << scoresFile << " does not exist." << endl;
		exit(1);
	}
}

void Coordinator::extractRepeatIndex(ScoresMaker * scoresMaker, int j) {
	vector<unsigned int> * rep = repeats->at(j);
	int * score = scoresMaker->getScores();
	int maxScore = scoresMaker->getMaxScore();
	int rowSize = 3 + maxScore + 1; // The one is for the zero score
	int str = segment->at(j)->at(0);

	bool inRpt = false;
	bool canFill = false;
	int s = -1;
	int e = -1;

	for (int i = 0; i < rep->size(); i++) {
		// Start a new MS
		if (rep->at(i) == 1 && !inRpt) {
			inRpt = true;
			s = i;
		}
		// End a the current MS
		else if (rep->at(i) == 0 && inRpt) {
			e = i - 1;
			inRpt = false;
			canFill = true;
		}
		// If the current MS at the end of the segment
		else if (i == rep->size() - 1 && inRpt) {
			e = i;
			inRpt = false;
			canFill = true;
		}
		// Extract features of the just recognized MS
		if (canFill) {
			vector<double> * index = new vector<double> ;
			index->push_back(s + str);
			index->push_back(e + str);
			index->push_back(getMean(score, s, e));
			repIndex->push_back(index);

			s = -1;
			e = -1;
			canFill = false;
		}
	}
}

void Coordinator::filter(vector<double> * paramList) {
	filter(paramList, 0.5);
}

void Coordinator::filter(vector<double> * paramList, double threshold) {
	if (!canFilter) {
		throw InvalidOrderOfOperationsException(
				"Filtering must occur after detection.");
	}

	double m1 = paramList->at(0);
	double m2 = paramList->at(1);
	double s1 = paramList->at(2);
	double s2 = paramList->at(3);
	double b = paramList->at(4);
	double w1 = paramList->at(5);
	double w2 = paramList->at(6);

	vector<vector<double> *> * copy = new vector<vector<double> *> ();

	for (int j = 0; j < repIndex->size(); j++) {
		double x1 = w1
				* ((repIndex->at(j)->at(1) - repIndex->at(j)->at(0) - m1) / s1);
		double x2 = w2 * ((repIndex->at(j)->at(2) - m2) / s2);
		double t = 1.0 / (1.0 + exp(-1.0 * (b + x1 + x2)));
		if (t >= threshold) {
			copy->push_back(new vector<double> (*(repIndex->at(j))));
		}
	}

	for (int j = 0; j < repIndex->size(); j++) {
		repIndex->at(j)->clear();
	}
	repIndex->clear();

	repIndex = copy;

	// Can output now
	canOutput = true;
}

void Coordinator::maskRepeats(string maskedFile) {
	if (!canOutput) {
		throw InvalidOrderOfOperationsException(
				"Writing must occur after detection or filtering.");
	}

	// Mask STR in sequence
	string baseCopy = base;
	for (int j = 0; j < repIndex->size(); j++) {
		for (int h = repIndex->at(j)->at(0); h <= repIndex->at(j)->at(1); h++) {
			baseCopy[h] = tolower(baseCopy[h]);
		}
	}

	// Write masked sequence in fasta format
	ofstream outMask;
	outMask.open(maskedFile.c_str(), ios::out | ios::app);

	if (!outMask) {
		cerr << "Error: " << maskedFile << " does not exist." << endl;
		exit(0);
	}

	outMask << info << endl;
	int step = 50;
	int len = baseCopy.size();
	for (int i = 0; i < len; i = i + step) {
		int e = (i + step - 1 > len - 1) ? len - 1 : i + step - 1;
		for (int k = i; k <= e; k++) {
			outMask << baseCopy[k];
		}
		outMask << endl;
	}

	outMask.close();
}

void Coordinator::writeRepeatsCoordinates(string indexFile) {
	if (!canOutput) {
		throw InvalidOrderOfOperationsException(
				"Writing must occur after detection or filtering.");
	}

	// Write repeat index
	ofstream outIndex;
	outIndex.open(indexFile.c_str(), ios::out | ios::app);

	// Make sure that the output file exists
	if (!outIndex) {
		cout << "Error: " << indexFile << " does not exist." << endl;
		exit(0);
	}

	// Write the index of the repeat segment [x,y[ with respect with the start (chrK:start-end)
	for (int j = 0; j < repIndex->size(); j++) {
		outIndex << chrm << ":" << ((int) (coorStart + repIndex->at(j)->at(0)))
				<< "-" << ((int) (coorStart + repIndex->at(j)->at(1) + 1))
				<< " ";
		for (int o = 2; o < repIndex->at(j)->size(); o++) {
			outIndex << repIndex->at(j)->at(o) << " ";
		}
		outIndex << endl;
	}
	// outIndex << endl;
	outIndex.close();
}

// Can be combined with the next method to speed the program
double Coordinator::getMean(int * scores, int start, int end) {
	double sum = 0;
	for (int i = start; i <= end; i++) {
		sum += scores[i];
	}
	double mean = sum / (double) (end - start + 1);
	return mean;
}

// Can be combined with the next method to speed the program
double Coordinator::getSum(int * scores, int start, int end) {
	double sum = 0;
	for (int i = start; i <= end; i++) {
		sum += scores[i];
	}

	return sum;
}
