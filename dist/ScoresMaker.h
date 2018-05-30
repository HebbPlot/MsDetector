/*
 * ScoresMaker.h
 *
 *  Created on: Oct 18, 2011
 *      Author: Hani Zakaria Girgis, Ph.D. at the NCBI/NLM/NIH - USA
 */

#ifndef SCORESMAKER_H_
#define SCORESMAKER_H_

#include <string>

using namespace std;

class ScoresMaker {
public:
	ScoresMaker(const char *, int, int, int, int);
	virtual ~ScoresMaker();
	void start();
	int * getScores();
	int getScoresSize();
	int getStart();
	int getMaxScore();

protected:
	const char * seq;
	int * scores;
	int s;
	int e;
	int wordLength;
	int w;
	int scoresSize;
	int maxScore;

	void check(int);
	virtual void score();
	virtual void scoreOppositDirection();
	virtual void fixPerfectMsBoundary();
};

#endif /* SCORESMAKER_H_ */
