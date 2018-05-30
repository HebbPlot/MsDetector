/*
 * TransitionCompositionMatrix.h
 *
 *  Created on: May 1, 2012
 *      Author: Hani Zakaria Girgis, Ph.D. at the NCBI/NLM/NIH - USA
 */

#ifndef TRANSITIONCOMPOSITIONMATRIX_H_
#define TRANSITIONCOMPOSITIONMATRIX_H_

#include "CompositionMatrix.h"

class TransitionCompositionMatrix: public CompositionMatrix {
public:
	TransitionCompositionMatrix(string);
	virtual ~TransitionCompositionMatrix();
protected:
	virtual void buildMatchTable();
};

#endif /* TRANSITIONCOMPOSITIONMATRIX_H_ */
