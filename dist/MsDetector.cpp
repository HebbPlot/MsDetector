//============================================================================
// Name        : MsDetector.cpp
// Author      : Hani Zakaria Girgis, Ph.D. at the NCBI/NLM/NIH - USA
//============================================================================

#include "Coordinator.h"
#include "Util.h"
#include "CompositionMatrix.h"
#include "TransitionCompositionMatrix.h"

#include <HMMlib/hmm_table.hpp>
#include <HMMlib/hmm_vector.hpp>
#include <HMMlib/hmm.hpp>

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <map>

using namespace std;

const static string SEQ_PRM = string("-seq");
const static string SCR_PRM = string("-scr");
const static string MSK_PRM = string("-msk");
const static string RPT_PRM = string("-rpt");
const static string HMM_PRM = string("-hmm");
const static string GLM_PRM = string("-glm");
const static string FRQ_PRM = string("-frq");
const static string FCT_PRM = string("-fct");
const static string MTR_PRM = string("-mtr");
const static string LEN_PRM = string("-len");
const static string SFL_PRM = string("-sfl");
const static string THR_PRM = string("-thr");
static CompositionMatrix * compMtr;
static TransitionCompositionMatrix * transCompMtr;

HMM<double> * buildHmm() {
	int K = 2; // number of states
	int M = 7; // alphabet size

	boost::shared_ptr<HMMVector<double> > pi_ptr(new HMMVector<double> (K));
	boost::shared_ptr<HMMMatrix<double> > T_ptr(new HMMMatrix<double> (K, K));
	boost::shared_ptr<HMMMatrix<double> > E_ptr(new HMMMatrix<double> (M, K));

	// initial probabilities
	HMMVector<double> &pi = *pi_ptr;
	double pSum = 39333 + 337;
	pi(0) = 39333 / pSum;
	pi(1) = 337 / pSum;
	// cout << pi(0) << " " << pi(1) << endl;

	// transitions
	HMMMatrix<double> &T = *T_ptr;
	double sumS0 = 19640711 + 2707;
	double sumS1 = 2657 + 149255;
	T(0, 0) = 19640711 / sumS0;
	T(0, 1) = 2707 / sumS0;
	T(1, 0) = 2657 / sumS1;
	T(1, 1) = 149255 / sumS1;

	// cout << T(0, 0) << " " << T(0, 1) << endl;
	// cout << T(1, 0) << " " << T(1, 1) << endl;

	HMMMatrix<double> &E = *E_ptr;
	double sumO0 = 0 + 53 + 35314 + 4041889 + 11972469 + 3398256 + 218675;
	E(0, 0) = 0 / sumO0;
	E(1, 0) = 53 / sumO0;
	E(2, 0) = 35314 / sumO0;
	E(3, 0) = 4041889 / sumO0;
	E(4, 0) = 11972469 / sumO0;
	E(5, 0) = 3398256 / sumO0;
	E(6, 0) = 218675 / sumO0;

	// Need to update the new numbers
	double sumO1 = 0 + 3 + 74 + 3246 + 19618 + 48184 + 97219;
	E(0, 1) = 0 / sumO1;
	E(1, 1) = 3 / sumO1;
	E(2, 1) = 74 / sumO1;
	E(3, 1) = 3246 / sumO1;
	E(4, 1) = 19618 / sumO1;
	E(5, 1) = 48184 / sumO1;
	E(6, 1) = 97219 / sumO1;

	// cout << E(0, 0) << " " << E(1, 0) << " " << E(2, 0) << " " << E(3, 0)
	//		<< " " << E(4, 0) << " " << E(5, 0) << " " << E(6, 0) << endl;

	// cout << E(0, 1) << " " << E(1, 1) << " " << E(2, 1) << " " << E(3, 1)
	//		<< " " << E(4, 1) << " " << E(5, 1) << " " << E(6, 1) << endl;

	// cout << "Constructing HMM" << endl;

	HMM<double> * myHMM = new HMM<double> (pi_ptr, T_ptr, E_ptr);
	return myHMM;
}

HMM<double> * buildHmmFromFile(string probFile) {
	ifstream in(probFile.c_str());

	int K = 2; // number of states

	boost::shared_ptr<HMMVector<double> > pi_ptr(new HMMVector<double> (K));
	boost::shared_ptr<HMMMatrix<double> > T_ptr(new HMMMatrix<double> (K, K));

	// initial probabilities
	HMMVector<double> &pi = *pi_ptr;
	double pSum = 0;

	for (int i = 0; i < K; i++) {
		in >> pi(i);
		pSum += pi(i);
	}

	for (int i = 0; i < K; i++) {
		pi(i) = pi(i) / pSum;
	}

	// cout << pi(0) << " " << pi(1) << endl;

	// transitions
	HMMMatrix<double> &T = *T_ptr;
	for (int i = 0; i < K; i++) {
		double sum = 0;
		for (int j = 0; j < K; j++) {
			in >> T(i, j);
			sum += T(i, j);
		}
		for (int j = 0; j < K; j++) {
			T(i, j) = T(i, j) / sum;
		}
	}

	// cout << T(0, 0) << " " << T(0, 1) << endl;
	// cout << T(1, 0) << " " << T(1, 1) << endl;

	// emission
	vector<double> * emissionTemp = new vector<double> ();
	double temp = -100;
	int counter = 1;
	while (in >> temp) {
		emissionTemp->push_back(temp);
	}
	in.close();

	int M = emissionTemp->size() / 2;
	// cout << "The number of outputs is: " << M << endl;

	boost::shared_ptr<HMMMatrix<double> > E_ptr(new HMMMatrix<double> (M, K));
	HMMMatrix<double> &E = *E_ptr;
	for (int i = 0; i < K; i++) {
		double sum = 0;
		for (int j = 0; j < M; j++) {
			E(j, i) = emissionTemp->at(i * M + j);
			sum += E(j, i);
		}
		for (int j = 0; j < M; j++) {
			E(j, i) = E(j, i) / sum;
		}
	}

	/*
	 for (int i = 0; i < K; i++) {
	 for (int j = 0; j < M; j++) {
	 cout << E(j, i) << " ";
	 }
	 cout << endl;
	 }
	 */

	emissionTemp->clear();
	delete emissionTemp;

	HMM<double> * myHMM = new HMM<double> (pi_ptr, T_ptr, E_ptr);

	return myHMM;
}

vector<double> * filter() {
	// Balanced
	/*
	 double m1 = 31.8856;
	 double m2 = 5.4510;
	 double s1 = 44.7814;
	 double s2 = 0.2939;
	 double b = 3.9949;
	 double w1 = 12.9101;
	 double w2 = 2.5953;
	 */

	// Unbalanced old
	/*
	 double m1 = 19.2408;
	 double m2 = 5.3453;
	 double s1 = 29.4133;
	 double s2 = 0.2682;
	 double b = -2.0210;
	 double w1 = 7.9237;
	 double w2 = 2.4688;
	 */

	// Current in the ISMB Paper
	/*
	 double m1 = 23.3066;
	 double m2 = 5.3375;
	 double s1 = 32.8712;
	 double s2 = 0.2712;
	 double b  = -1.3475;
	 double w1 = 8.0199;
	 double w2 = 2.4840;
	 */

	// Similar Current in the ISMB Paper (You should use this!)
	/*
	 double m1 = 22.0027;
	 double m2 = 5.3350;
	 double s1 = 27.7575;
	 double s2 = 0.2692;
	 double b = -1.7158;
	 double w1 = 6.7572;
	 double w2 = 2.4587;
	 */

	/*
	 double m1 = 21.0727;
	 double m2 = 5.2383;
	 double s1 = 20.6329;
	 double s2 = 0.2853;
	 double b  = -2.1134;
	 double w1 = 5.2373;
	 double w2 = 2.6218;
	 */

	/*
	 double m1 = 23.5779;
	 double m2 = 5.4002;
	 double s1 = 23.9167;
	 double s2 = 0.2073;
	 double b  = -0.9538;
	 double w1 = 5.7861;
	 double w2 = 2.4425;
	 */

	// 10K
	/*
	 double m1 = 33.1880;
	 double m2 = 5.4230;
	 double s1 = 39.1042;
	 double s2 = 0.3054;
	 double b  = 3.5541;
	 double w1 = 11.6120;
	 double w2 = 3.2277;
	 */

	// 1K
	/*
	 double m1 = 22.3743;
	 double m2 = 5.3384;
	 double s1 = 28.2782;
	 double s2 = 0.2702;
	 double b  = -1.5576;
	 double w1 = 7.5018;
	 double w2 = 2.6757;
	 */

	double m1 = 20.909;
	double m2 = 5.327;
	double s1 = 26.208;
	double s2 = 0.26344;
	double b = -2.2456;
	double w1 = 8.2561;
	double w2 = 2.9729;

	vector<double> * paramList = new vector<double> ();
	paramList->push_back(m1);
	paramList->push_back(m2);
	paramList->push_back(s1);
	paramList->push_back(s2);
	paramList->push_back(b);
	paramList->push_back(w1);
	paramList->push_back(w2);

	return paramList;
}

vector<double> * buildFilterFromFile(string filterFile) {
	ifstream in(filterFile.c_str());

	vector<double> * paramList = new vector<double> ();
	double temp;
	while (in >> temp) {
		paramList->push_back(temp);
	}
	in.close();

	return paramList;
}

void coordinate(map<string, string>* param) {
	// Read sequence file
	vector<string> * infoList = new vector<string> ();
	infoList->reserve(1000);
	vector<string> * seqList = new vector<string> ();
	seqList->reserve(1000);
	Util::readFasta(param->at(SEQ_PRM), infoList, seqList);

	// Construct HMM
	// HMM<double> * hmmPtr = buildHmm();
	HMM<double> * hmmPtr;
	if (param->count(HMM_PRM) > 0) {
		Util::checkFile(param->at(HMM_PRM));
		hmmPtr = buildHmmFromFile(param->at(HMM_PRM));
	}

	// Construct Filter
	vector<double> * filterPtr;
	if (param->count(GLM_PRM) > 0) {
		Util::checkFile(param->at(GLM_PRM));
		filterPtr = buildFilterFromFile(param->at(GLM_PRM));
	}

	bool canShuffle = false;
	if (param->count(SFL_PRM) > 0) {
		string value = param->at(SFL_PRM);
		if (value.compare("1") == 0) {
			canShuffle = true;
		}
	}

	if (param->count(FRQ_PRM) > 0) {
		Util::checkFile(param->at(FRQ_PRM));

		compMtr = new CompositionMatrix(param->at(FRQ_PRM));
		compMtr->start();
		// compMtr->printTable();

		transCompMtr = new TransitionCompositionMatrix(param->at(FRQ_PRM));
		transCompMtr->start();
		// transCompMtr->printTable();
	}

	// Delete output files if they exist.
	if (param->count(SCR_PRM) > 0) {
		Util::deleteFile(param->at(SCR_PRM));
	}
	if (param->count(MSK_PRM) > 0) {
		Util::deleteFile(param->at(MSK_PRM));

	}
	if (param->count(RPT_PRM) > 0) {
		Util::deleteFile(param->at(RPT_PRM));
	}

	for (int i = 0; i < infoList->size(); i++) {
		Coordinator * coor = new Coordinator(seqList->at(i), infoList->at(i),
				atoi((param->at(LEN_PRM)).c_str()), atoi(
						(param->at(FCT_PRM)).c_str()), canShuffle);

		// Operations
		// 1-a: scoring
		string mtrx = param->at(MTR_PRM);
		if (mtrx.compare(string("Comp")) == 0) {
			coor->score(compMtr);
		} else if (mtrx.compare(string("TransComp")) == 0) {
			coor->score(transCompMtr);
		} else {
			coor->score(mtrx);
		}

		// 1-b: scores writing
		if (param->count(SCR_PRM) > 0) {
			coor->writeScores(param->at(SCR_PRM));
		}
		// 2: detection
		if (param->count(HMM_PRM) > 0) {
			coor->detect(hmmPtr);
		}
		// 3: filtering
		if (param->count(GLM_PRM) > 0) {
			if (param->count(THR_PRM) > 0) {
				coor->filter(filterPtr, atof(param->at(THR_PRM).c_str()));
			} else {
				coor->filter(filterPtr);
			}
		}
		// 4-a: masking repeats
		if (param->count(MSK_PRM) > 0) {
			coor->maskRepeats(param->at(MSK_PRM));
		}
		// 4-b: writing repeats coordinates
		if (param->count(RPT_PRM) > 0) {
			coor->writeRepeatsCoordinates(param->at(RPT_PRM));
		}

		delete coor;
	}

	// Free memory
	if (param->count(HMM_PRM) > 0) {
		delete hmmPtr;
	}

	if (param->count(GLM_PRM) > 0) {
		filterPtr->clear();
		delete filterPtr;
	}

	infoList->clear();
	delete infoList;
	seqList->clear();
	delete seqList;
}

void printMassage(){
        cout << endl;
        cout << "MsDetector is designed and developed by Hani Zakaria Girgis under the supervision" << endl;
        cout << "of Sergey Sheetlin at the Spouge Research Group, the National Center for" << endl;
        cout << "Biotechnology Information, The National Institutes of Health, USA." << endl;
        cout << endl;
        cout << "The library HMMlib is utilized in MsDetector. HMMlib is covered by the GNU Lesser" << endl;
        cout << "General Public License." << endl;
        cout << endl;
        cout << "Version: 1.2" << endl;
        cout << endl;
}

int main(int argc, char * argv[]) {
	printMassage();
	// Usage massage
	string massage = string("Valid arguments pairs:\n");
	massage.append("\t-seq sequence_file\n");
	massage.append("\t-sfl shuffle [0 or 1 default is 0]\n");
	massage.append("\t-hmm hmm_file\n");
	massage.append("\t-glm glm_file\n");
	massage.append("\t-thr glm_threshold [>= 0.0 && <= 1.0]");
	massage.append("\t-len word_length\n");
	massage.append("\t-fct win_factor [window = 2 x fct x len]\n");
	massage.append("\t-mtr matrix_name [Id, Trans, Comp, and TransComp]\n");
	massage.append(
			"\t-frq frequencies_file [must be provided with the Comp or the TransComp matrixes]\n");
	massage.append("\t-scr scores_file\n");
	massage.append("\t-msk masked_sequence_file\n");
	massage.append("\t-rpt repeats_index_file\n");

	// Table of valid argument pairs
	map<string, string> * validParam = new map<string, string> ();
	validParam->insert(map<string, string>::value_type(SEQ_PRM, "DUMMY"));
	validParam->insert(map<string, string>::value_type(SFL_PRM, "DUMMY"));
	validParam->insert(map<string, string>::value_type(SCR_PRM, "DUMMY"));
	validParam->insert(map<string, string>::value_type(MSK_PRM, "DUMMY"));
	validParam->insert(map<string, string>::value_type(RPT_PRM, "DUMMY"));
	validParam->insert(map<string, string>::value_type(HMM_PRM, "DUMMY"));
	validParam->insert(map<string, string>::value_type(GLM_PRM, "DUMMY"));
	validParam->insert(map<string, string>::value_type(FRQ_PRM, "DUMMY"));
	validParam->insert(map<string, string>::value_type(FCT_PRM, "DUMMY"));
	validParam->insert(map<string, string>::value_type(MTR_PRM, "DUMMY"));
	validParam->insert(map<string, string>::value_type(LEN_PRM, "DUMMY"));
	validParam->insert(map<string, string>::value_type(THR_PRM, "DUMMY"));

	// Make a table of the user provided arguments
	map<string, string> * param = new map<string, string> ();
	if (argc > 1 && argc % 2 == 1) {
		for (int i = 1; i < argc - 1; i += 2) {
			if (validParam->count(argv[i]) > 0) {
				param->insert(map<string, string>::value_type(argv[i], argv[i
						+ 1]));
			} else {
				cerr << "Invalid argument: " << argv[i] << " " << argv[i + 1]
						<< endl;
				cerr << massage << endl;
				return 0;
			}
		}

		// Check if the user provided the essential arguments
		if (param->count(SEQ_PRM) == 0) {
			cerr << "Please provide the input sequence file." << endl;
			cerr << massage << endl;
			return 0;
		} else {
			Util::checkFile(param->at(SEQ_PRM));
		}

		if (param->count(LEN_PRM) == 0) {
			cerr << "Please provide the word length." << endl;
			cerr << massage << endl;
			return 0;
		}

		if (param->count(FCT_PRM) == 0) {
			cerr << "Please provide the window factor." << endl;
			cerr << massage << endl;
			return 0;
		}

		if (param->count(MTR_PRM) == 0) {
			cerr << "Please provide a matrix name." << endl;
			cerr << massage << endl;
			return 0;
		}

		if (param->count(MTR_PRM) > 0) {
			string mtrx = param->at(MTR_PRM);
			if (mtrx.compare(string("Comp")) == 0 || mtrx.compare(string(
					"TransComp")) == 0) {
				if (param->count(FRQ_PRM) == 0) {
					cerr << "Please provide the frequency file." << endl;
					cerr << massage << endl;
					return 0;
				}
			} else if (mtrx.compare(string("Id")) == 0 || mtrx.compare(string(
					"Trans")) == 0) {
				// it is a valid matrix
			} else {
				cerr << "Please provide a valid matrix name." << endl;
				cerr << massage << endl;
				return 0;
			}
		}

		coordinate(param);
	} else {
		cerr << "Please provide argument pairs of the form: -flag value."
				<< endl;
		cerr << massage << endl;
	}

	// Clean up
	validParam->clear();
	delete validParam;
	param->clear();
	delete param;

	return 0;
}
