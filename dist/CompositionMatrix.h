/*
 * ProbabilitiesReader.h
 *
 *  Created on: May 1, 2012
 *      Author: Hani Zakaria Girgis, Ph.D. at the NCBI/NLM/NIH - USA
 */

#ifndef COMPOSITIONMATRIX_H_
#define COMPOSITIONMATRIX_H_

#include <string>
#include <vector>
#include <map>

using namespace std;

class CompositionMatrix {
public:
	static const int indexShift;

	CompositionMatrix(string);
	virtual ~CompositionMatrix();
	map<char, double> * getProbTable();
	map<int, double> * getMatchTable();
	double getBestMatchScore();
	void printTable();
	void start();
	double * getMatchList();

protected:
	string probFile;
	static const int matchListSize;
	map<char, double> * probTable;
	map<int, double> * matchTable;
	double bestMatchScore;
	vector<char> * bases;
	double * matchList;

	void buildProbTable();
	virtual void buildMatchTable();
	void findBestMatch();
	void buildMatchList();
};

#endif /* PROBABILITIESREADER_H_ */
