#ifndef ALGORITHMCOMPARISON_H
#define ALGORITHMCOMPARISON_H

#include "cvrlrocStructs.h"
#include "ROC.h"
#include "ScoreData.h"

class AlgorithmComparison {
public:
	//data members holding the scores for the two algorithms
	ScoreData* a1;
	ScoreData* a2;
	//the associated ROCs with the algorithms
	ROC* roc1;
	ROC* roc2;
	//store the mean and stddev for the stats
	double a1_TAmean, a1_TAstd, a2_TAmean, a2_TAstd;
	//the correlation coefficient
	double cc;
	//the Z statistic
	double zStat;
	//the set FA rate
	double FA_target_rate;

	//default constructor
	AlgorithmComparison();
	//constructor that accepts score data and calculates everything
	AlgorithmComparison(ScoreData& alg1, ScoreData& alg2, int numBoot, int numROCPoint, double FApoint);

	//create bootstraps of corresponding scores
	bool createSynchronizedBootstraps(int num);	

	//calculate the correlateion for TA @ FA
	void calculateCorrelationTAatFA(double FApoint);
	//calculate the Z statistic between the two algorithms for TA @ FA
	void calculateZStatisticTAatFA();

	void printStatsToR(std::string filename,std::string label1, std::string label2);
};

#endif
