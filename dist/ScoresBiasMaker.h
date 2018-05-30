/*
 * ScoresBiasMaker.h
 *
 *  Created on: Apr 17, 2012
 *      Author: Hani Zakaria Girgis, Ph.D. at the NCBI/NLM/NIH - USA
 */

#ifndef SCORESBIASMAKER_H_
#define SCORESBIASMAKER_H_

#include "ScoresMaker.h"
#include "CompositionMatrix.h"

#include <map>
#include <math.h>

class ScoresBiasMaker: public ScoresMaker {
public:
	ScoresBiasMaker(const char *, int, int, int, int, CompositionMatrix *);
	virtual ~ScoresBiasMaker();

protected:
	map<char, double> * probTable;
	map<int, double> * matchTable;
	double * matchList;

	virtual void score();

	/**
	 * Credit: http://stackoverflow.com/questions/554204/where-is-round-in-c
	 */
	inline double round(double number) {
		return number < 0.0 ? ceil(number - 0.5) : floor(number + 0.5);
	}

private:
	/**
	 * The match between abundant nucleotide is less than that of less abundant nucleotide
	 */
	inline double probMatch(char b1, char b2) {
		double result = 0;

		if (probTable->count(b1) > 0 && probTable->count(b2) > 0) {
			// result = matchTable->at(b1 + b2);
			result = matchList[b1 + b2];
		}

		return result;
	}

};

#endif /* SCORESBIASMAKER_H_ */
