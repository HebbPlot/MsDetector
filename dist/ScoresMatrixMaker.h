/*
 * ScoresMatrixMaker.h
 *
 *  Created on: Apr 16, 2012
 *      Author: Hani Zakaria Girgis, Ph.D. at the NCBI/NLM/NIH  - USA
 */

#ifndef SCORESMATRIXMAKER_H_
#define SCORESMATRIXMAKER_H_

#include "ScoresMaker.h"

class ScoresMatrixMaker: public ScoresMaker {
public:
	ScoresMatrixMaker(const char *, int, int, int, int);
	virtual ~ScoresMatrixMaker();

private:
	virtual void score();

	/**
	 * Penalize transitions (1-ring <--> 1-ring and 2-ring <--> 2-ring)
	 * less than transversions (1-ring <--> 2-ring).
	 * Transitions include:
	 * 	A <--> G
	 * 	C <--> T
	 *
	 * b1: the first base
	 * b2: the second base
	 */
	inline double transMatch(char b1, char b2) {
		double result = 0.0;

		if (b1 == b2) {
			result = 2.0;
		} else {
			switch (b1) {
			case 'A': {
				if (b2 == 'G') {
					result = 1.0;
				}
			}
				break;
			case 'G': {
				if (b2 == 'A') {
					result = 1.0;
				}
			}
				break;
			case 'C': {
				if (b2 == 'T') {
					result = 1.0;
				}
			}
				break;
			case 'T': {
				if (b2 == 'C') {
					result = 1.0;
				}
			}
				break;
			default: {
				result = 0.0;
			}
				break;
			}
		}
		return result;
	}
};

#endif /* SCORESMATRIXMAKER_H_ */
