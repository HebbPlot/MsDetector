/*
 * ProbabilitiesReader.cpp
 *
 *  Created on: May 1, 2012
 *      Author: Hani Zakaria Girgis, Ph.D. at the NCBI/NLM/NIH - USA
 */

#include "CompositionMatrix.h"
#include "Util.h"
#include "InvalidInputException.h"

#include <string>
#include <vector>
#include <map>
#include <math.h>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

const int CompositionMatrix::matchListSize = 1 + 'T' + 'T';

CompositionMatrix::CompositionMatrix(string bFile) {
	probFile = bFile;

	bases = new vector<char> ();
	bases->push_back('A');
	bases->push_back('C');
	bases->push_back('G');
	bases->push_back('T');

	probTable = new map<char, double> ();
	matchTable = new map<int, double> ();

	matchList = new double[matchListSize];
	for (int i = 0; i < matchListSize; i++) {
		matchList[i] = 0.0;
	}
}

CompositionMatrix::~CompositionMatrix() {
	probTable->clear();
	delete probTable;

	matchTable->clear();
	delete matchTable;

	delete matchList;
}

void CompositionMatrix::start() {
	buildProbTable();
	buildMatchTable();
	findBestMatch();
	buildMatchList();
}

void CompositionMatrix::buildProbTable() {
	// Read probabilities from file
	ifstream in(probFile.c_str());
	double a;
	double c;
	double g;
	double t;
	in >> a;
	in >> c;
	in >> g;
	in >> t;
	in.close();

	// Save probabilities into a hash table
	probTable->insert(map<char, double>::value_type('A', a));
	probTable->insert(map<char, double>::value_type('C', c));
	probTable->insert(map<char, double>::value_type('G', g));
	probTable->insert(map<char, double>::value_type('T', t));

	// Check these base probabilities
	map<char, double>::iterator iter = probTable->begin();
	map<char, double>::iterator endIter = probTable->end();

	double sum = 0.0;
	while (iter != endIter) {
		char key = iter->first;
		double value = iter->second;
		sum += value;
		if (value < 0 || value > 1) {
			string msg("Invalid probabilities found in file: ");
			msg.append(probFile);
			msg.append("\n");
			msg.append("Base ");
			const char keyArr[] = { key };
			msg.append(keyArr);
			msg.append(" has an invalid probability of ");
			msg.append(Util::double2string(value));
			throw InvalidInputException(msg);
		}
		iter++;
	}

	if (abs(1.0 - sum) > 0.000001) {
		string msg("Invalid sum of probabilities found in file: ");
		msg.append(probFile);
		msg.append("\n");
		msg.append("Base probabilities must add up to 1. The current sum is ");
		msg.append(Util::double2string(sum));
		throw InvalidInputException(msg);
	}
}
void CompositionMatrix::buildMatchTable() {
	// Calculate the score corrected according to the nucleotide composition
	// The match score is positive and the mismatch score is negative.
	for (int i = 0; i < bases->size(); i++) {
		for (int j = i; j < bases->size(); j++) {
			int key = bases->at(i) + bases->at(j);
			double value = 0.0;
			if (bases->at(i) == bases->at(j)) {
				value = -1.0 * log(probTable->at(bases->at(i)));
			} else {
				value = 0.5 * (log(probTable->at(bases->at(i))) + log(
						probTable->at(bases->at(j))));
			}
			matchTable->insert(map<int, double>::value_type(key, value));
		}
	}
}

void CompositionMatrix::findBestMatch() {
	// Find worst match
	bestMatchScore = -1000000;
	map<int, double>::iterator iter = matchTable->begin();
	map<int, double>::iterator iterEnd = matchTable->end();
	while (iter != iterEnd) {
		if (iter->second > bestMatchScore) {
			bestMatchScore = iter->second;
		}
		iter++;
	}
}

void CompositionMatrix::buildMatchList() {
	map<int, double>::iterator iter = matchTable->begin();
	map<int, double>::iterator iterEnd = matchTable->end();
	while (iter != iterEnd) {
		matchList[iter->first] = iter->second;
		iter++;
	}
}

void CompositionMatrix::printTable() {
	map<int, double>::iterator iter = matchTable->begin();
	map<int, double>::iterator iterEnd = matchTable->end();
	while (iter != iterEnd) {
		cout << iter->first << " -> " << iter->second << endl;
		iter++;
	}
}

double CompositionMatrix::getBestMatchScore() {
	return bestMatchScore;
}

map<char, double> * CompositionMatrix::getProbTable() {
	return probTable;
}

map<int, double> * CompositionMatrix::getMatchTable() {
	return matchTable;
}

double * CompositionMatrix::getMatchList() {
	return matchList;
}
