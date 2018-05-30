/*
 * ScoresMatrixMaker.cpp
 *
 *  Created on: Apr 16, 2012
 *      Author: Hani Zakaria Girgis, Ph.D. at the NCBI/NLM/NIH - USA
 */
#include "ScoresMatrixMaker.h"
#include <iostream>

ScoresMatrixMaker::ScoresMatrixMaker(const char * sequence, int start, int end,
		int l, int winLength) :
	ScoresMaker(sequence, start, end, l, winLength) {
	maxScore = 2 * wordLength;
}

ScoresMatrixMaker::~ScoresMatrixMaker() {
	// TODO Auto-generated destructor stub
}

void ScoresMatrixMaker::score() {
	// First n words. n is the length of the word
	for (int i = s; i < s + wordLength; i++) {
		int score1 = -1;
		int rStart = i + 1;
		int rEnd = rStart + w - wordLength + 1;
		rEnd = (rEnd > e - wordLength + 1) ? e - wordLength + 1 : rEnd;

		for (int r = rStart; r <= rEnd; r++) {
			double sum = 0;
			for (int j = 0; j < wordLength; j++) {
				sum += transMatch(seq[i + j], seq[r + j]);
			}
			score1 = (sum > score1) ? sum : score1;
			if (score1 == maxScore) {
				break;
			}
		}
		scores[i - s] = score1;
	}

	// Main body
	for (int i = s + wordLength; i <= e - wordLength; i++) {
		// Left
		int lStart = (i - w > s) ? i - w : s;
		int lEnd = i - 1;
		int score2 = -1;
		for (int l = lStart; l <= lEnd; l++) {
			double sum = 0.0;
			for (int j = 0; j < wordLength; j++) {
				sum += transMatch(seq[i + j], seq[l + j]);
			}
			score2 = (sum > score2) ? sum : score2;
			if (score2 == maxScore) {
				break;
			}
		}

		// Right
		if (score2 < maxScore) {
			int rStart = i + 1;
			int rEnd = rStart + w - wordLength + 1;
			rEnd = (rEnd > e - wordLength + 1) ? e - wordLength + 1 : rEnd;

			for (int r = rStart; r <= rEnd; r++) {
				double sum = 0.0;
				for (int j = 0; j < wordLength; j++) {
					sum += transMatch(seq[i + j], seq[r + j]);
				}
				score2 = (sum > score2) ? sum : score2;
				if (score2 == maxScore) {
					break;
				}
			}
		}

		scores[i - s] = score2;
	}

	// Last word
	int i = e - wordLength + 1;
	int lStart = (i - w > s) ? i - w : s;
	int lEnd = i - 1;
	int score3 = -1;
	for (int l = lStart; l <= lEnd; l++) {
		double sum = 0.0;
		for (int j = 0; j < wordLength; j++) {
			sum += transMatch(seq[i + j], seq[l + j]);
		}
		score3 = (sum > score3) ? sum : score3;
		if (score3 == maxScore) {
			break;
		}
	}
	scores[i - s] = score3;

	for (int h = e - wordLength + 2 - s; h <= e - s; h++) {
		scores[h] = scores[i - s];
	}
}
