/*
 * TransitionCompositionMatrix.cpp
 *
 *  Created on: May 1, 2012
 *      Author: Hani Zakaria Girgis, Ph.D. at the NCBI/NLM/NIH - USA
 */

#include "TransitionCompositionMatrix.h"
#include "CompositionMatrix.h"

#include <math.h>
#include <iostream>

using namespace std;

TransitionCompositionMatrix::TransitionCompositionMatrix(string pFile) :
	CompositionMatrix(pFile) {
}

TransitionCompositionMatrix::~TransitionCompositionMatrix() {

}

void TransitionCompositionMatrix::buildMatchTable() {
	CompositionMatrix::buildMatchTable();
	(*matchTable)['A' + 'G'] = -0.25 * (log(probTable->at('A')) + log(
			probTable->at('G')));
	(*matchTable)['C' + 'T'] = -0.25 * (log(probTable->at('C')) + log(
			probTable->at('T')));
}
