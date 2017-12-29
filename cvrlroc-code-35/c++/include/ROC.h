#ifndef ROC_H
#define ROC_H

#include "cvrlrocStructs.h"
#include "ScoreData.h"
#include "Parameters.h"

class ROC {
public:
	//declare vector to hold TA,FA at the various points
	std::vector<std::vector<ROCpoint> > points;
	//declare vector hold the mean and std deviation of each point
	std::vector<ROCstat> stats;
	//declare a vector of confidence intervals
	std::vector<CI> CIs;
	//pointer to underlying ScoreData object
	ScoreData* scoredata;
	//should be a sorted threshold container
	std::vector<double> thresholds;
	int numBoot;

	//Default constructor
	ROC();
	//constructor to use score data to compute ROCs with numBoot bootstraps
	ROC(ScoreData& scores, Parameters& p);

	//calculate ROCs for a vector of vector of scores
	void calculateROC_bootstrap(ScoreData& scores, Parameters& p);

	//calculate a single ROC point by point
	void calculateSingleROC(const std::vector<double>& matchScores, const std::vector<double>& nonmatchScores);

	//calculate stats from curves
	void calculateStats();
	//print the statistics
	void printStats();
	//print in R format to a file
	void printStatsToR(Parameters& p);

	//function to compare doubles
	//returns true if they are equal, false otherwise
	bool dbleq(double a, double b);
};

#endif
