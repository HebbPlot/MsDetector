/*
 * ScoresMaker.cpp
 *
 *  Created on: Oct 18, 2011
 *      Author: Hani Zakaria Girgis, Ph.D. at the NCBI/NLM/NIH - USA
 */
#include "ScoresMaker.h"
#include "InvalidScoreException.h"
#include "Util.h"

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

ScoresMaker::ScoresMaker(const char * sequence, int start, int end, int l,
		int winLength) {
	seq = sequence;
	s = start;
	e = end;
	wordLength = l;
	maxScore = wordLength;
	w = winLength;
	// scoresSize = e - s + 1 - wordLength + 1;
	// scores = new int[e - s + 1];

	scoresSize = e - s + 1;
	scores = new int[scoresSize];
	for (int i = 0; i < scoresSize; i++) {
		scores[i] = -1;
	}
}

void ScoresMaker::start() {
	score();
	// scoreOppositDirection();
	// This can be commented out
	fixPerfectMsBoundary();
	check(scoresSize);
}

ScoresMaker::~ScoresMaker() {
	delete[] scores;
}

void ScoresMaker::score() {
	// First n words. n is the length of the word
	for (int i = s; i < s + wordLength; i++) {
		int score1 = -1;
		//int rStart = i + wordLength;
		int rStart = i + 1;
		int rEnd = rStart + w - wordLength + 1;
		rEnd = (rEnd > e - wordLength + 1) ? e - wordLength + 1 : rEnd;

		for (int r = rStart; r <= rEnd; r++) {
			int sum = 0;
			for (int j = 0; j < wordLength; j++) {
				if (seq[i + j] == seq[r + j]) {
					sum++;
				}
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
		// int lEnd = i - wordLength;
		int lEnd = i - 1;
		int score2 = -1;
		for (int l = lStart; l <= lEnd; l++) {
			int sum = 0;
			for (int j = 0; j < wordLength; j++) {
				if (seq[i + j] == seq[l + j]) {
					sum++;
				}
			}
			score2 = (sum > score2) ? sum : score2;
			if (score2 == maxScore) {
				break;
			}
		}

		// Right
		if (score2 != maxScore) {
			// int rStart = i + wordLength;
			int rStart = i + 1;
			int rEnd = rStart + w - wordLength + 1;
			rEnd = (rEnd > e - wordLength + 1) ? e - wordLength + 1 : rEnd;

			for (int r = rStart; r <= rEnd; r++) {
				int sum = 0;
				for (int j = 0; j < wordLength; j++) {
					if (seq[i + j] == seq[r + j]) {
						sum++;
					}
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
	// int lEnd = i - wordLength;
	int lEnd = i - 1;
	int score3 = -1;
	for (int l = lStart; l <= lEnd; l++) {
		int sum = 0;
		for (int j = 0; j < wordLength; j++) {
			if (seq[i + j] == seq[l + j]) {
				sum++;
			}
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

// This method compares the scores at the start and the end of the word.
// If the score at the end is smaller than that at the start, the method
// swap them. The goal of the operation is to correct the boundaries of
// MSs.
void ScoresMaker::scoreOppositDirection() {
	int * tempScores = new int[scoresSize];
	for (int i = 0; i < scoresSize; i++) {
		tempScores[i] = scores[i];
	}

	for (int e = wordLength - 1; e < scoresSize; e++) {
		int s = e - wordLength + 1;
		if (tempScores[e] < tempScores[s]) {
			scores[e] = tempScores[s];
		}
	}

	delete[] tempScores;
}

void ScoresMaker::fixPerfectMsBoundary() {
	// Work on a copy of the scores
	int * tempScores = new int[scoresSize];
	for (int i = 0; i < scoresSize; i++) {
		tempScores[i] = scores[i];
	}

	int islandLength = 15;

	// Calculate the maximum score of the island
	// In this case, six consecutive maximum scores
	int islandScore = islandLength * maxScore;

	// 5 * 4
	int boundScore = 20;

	// Calculate the score of the current island
	int islandSum = 0;
	for (int j = 0; j < islandLength; j++) {
		islandSum += tempScores[j];
	}

	// Calculate the boundary sum
	int boundSum = 0;
	for (int b = islandLength; b < islandLength + wordLength - 1; b++) {
		boundSum += tempScores[b];
	}

	// Adjust boundary of perfect MSs
	for (int h = islandLength; h <= scoresSize - wordLength; h++) {
		// Found potential perfect MS boundary
		if (islandSum == islandScore && boundSum < boundScore) {
			// Testing code: print the island and the corrected boundary as well as the location
			// for (int y = h - islandLength; y < h + wordLength - 1; y++) {
			//	cout << seq[s+y];
			// }
			// cout << "\t";
			// cout << h - islandLength << "-" << h + wordLength - 1 << endl;
			// End testing code

			// Adjust boundary score
			for (int k = h; k < h + wordLength - 1; k++) {
				scores[k] = maxScore;
			}
		}

		// Update score; subtract the score at the start of the island and
		// add the score at the current index
		islandSum = islandSum - tempScores[h - islandLength] + tempScores[h];
		boundSum = boundSum - tempScores[h] + tempScores[h + wordLength - 1];
	}

	// Free memory
	delete[] tempScores;
}

int * ScoresMaker::getScores() {
	return scores;
}

void ScoresMaker::check(int size) {
	for (int j = 0; j < size; j++) {
		if (scores[j] < 0 || scores[j] > maxScore) {
			string massage("ScoresMaker: Wrong value of " + Util::int2string(
					scores[j]) + " at index: " + Util::int2string(j));
			throw InvalidScoreException(massage);
		}
	}
}

int ScoresMaker::getScoresSize() {
	return scoresSize;
}

int ScoresMaker::getStart() {
	return s;
}

int ScoresMaker::getMaxScore(){
	return maxScore;
}
