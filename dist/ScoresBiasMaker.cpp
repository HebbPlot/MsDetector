/*
 * ScoresBiasMaker.cpp
 *
 *  Created on: Apr 17, 2012
 *      Author: Hani Zakaria Girgis, Ph.D. at the NCBI/NLM/NIH - USA
 */
#include "ScoresBiasMaker.h"
#include "ScoresMaker.h"
#include "CompositionMatrix.h"
#include "InvalidScoreException.h"

#include <map>
#include <iostream>
#include <math.h>

using namespace std;

ScoresBiasMaker::ScoresBiasMaker(const char * sequence, int start, int end,
		int l, int winLength, CompositionMatrix * matrix) :
	ScoresMaker(sequence, start, end, l, winLength) {

	probTable = matrix->getProbTable();
	matchTable = matrix->getMatchTable();
	matchList = matrix->getMatchList();
	maxScore = round(wordLength * matrix->getBestMatchScore());
}

ScoresBiasMaker::~ScoresBiasMaker() {

}

/**
 * This method is a copy of the one found in the base class.
 * I need to do this to inline the proMatch method.
 * But, there should be a better way to do this.
 */
void ScoresBiasMaker::score() {
	// First n words. n is the length of the word
	for (int i = s; i < s + wordLength; i++) {
		double score1 = -1;
		// int rStart = i + wordLength;
		int rStart = i + 1;
		int rEnd = rStart + w - wordLength + 1;
		rEnd = (rEnd > e - wordLength + 1) ? e - wordLength + 1 : rEnd;

		for (int r = rStart; r <= rEnd; r++) {
			double sum = 0;
			for (int j = 0; j < wordLength; j++) {
				sum += probMatch(seq[i + j], seq[r + j]);
			}
			score1 = (sum > score1) ? sum : score1;
			if (((int) score1) == maxScore) {
				break;
			}
		}

		if (score1 > maxScore) {
			score1 = maxScore;
		} else if (score1 < 0.0) {
			score1 = 0.0;
		} else {
			score1 = round(score1);
		}
		scores[i - s] = score1;

	}

	// Main body
	for (int i = s + wordLength; i <= e - wordLength; i++) {
		// Left
		int lStart = (i - w > s) ? i - w : s;
		// int lEnd = i - wordLength;
		int lEnd = i - 1;
		double score2 = -1;
		for (int l = lStart; l <= lEnd; l++) {
			double sum = 0.0;
			for (int j = 0; j < wordLength; j++) {
				sum += probMatch(seq[i + j], seq[l + j]);
			}
			score2 = (sum > score2) ? sum : score2;
			if (((int) score2) == maxScore) {
				break;
			}
		}

		// Right
		if (score2 < maxScore) {
			// int rStart = i + wordLength;
			int rStart = i + 1;
			int rEnd = rStart + w - wordLength + 1;
			rEnd = (rEnd > e - wordLength + 1) ? e - wordLength + 1 : rEnd;

			for (int r = rStart; r <= rEnd; r++) {
				double sum = 0.0;
				for (int j = 0; j < wordLength; j++) {
					sum += probMatch(seq[i + j], seq[r + j]);
				}
				score2 = (sum > score2) ? sum : score2;
				if (((int) score2) == maxScore) {
					break;
				}
			}
		}

		if (score2 > maxScore) {
			score2 = maxScore;
		} else if (score2 < 0.0) {
			score2 = 0.0;
		} else {
			score2 = round(score2);
		}
		scores[i - s] = score2;
	}

	// Last word
	int i = e - wordLength + 1;
	int lStart = (i - w > s) ? i - w : s;
	// int lEnd = i - wordLength;
	int lEnd = i - 1;
	double score3 = -1;
	for (int l = lStart; l <= lEnd; l++) {
		double sum = 0.0;
		for (int j = 0; j < wordLength; j++) {
			sum += probMatch(seq[i + j], seq[l + j]);
		}
		score3 = (sum > score3) ? sum : score3;
		if (((int) score3) == maxScore) {
			break;
		}
	}

	if (score3 > maxScore) {
		score3 = maxScore;
	} else if (score3 < 0.0) {
		score3 = 0.0;
	} else {
		score3 = round(score3);
	}
	scores[i - s] = score3;

	for (int h = e - wordLength + 2 - s; h <= e - s; h++) {
		scores[h] = scores[i - s];
	}
}
