#ifndef SCOREDATA_H
#define SCOREDATA_H

#include "cvrlrocStructs.h"
#include "MersenneTwister.h"

class ScoreData {
public:
	//declare a vector<vector> of match scores and nonmatch scores
	std::vector<std::vector<double> > match;	
	std::vector<std::vector<double> > nonmatch;
	std::vector<double> matchDist;
	std::vector<double> nonmatchDist;
	std::vector<double> binMax;
	
	//how many bootstraps of the original data exist
	//this value is 1 if only the original data is present
	int numBoot;
	MTRand RNG;

	//default constructor to set everything to 0
	ScoreData();
	//constructor to load score data from file
	ScoreData(const std::string filename);
	//copy constructor
	ScoreData(const ScoreData& obj);

	//convert higher is better to lower is better
	void convertHigherToLower();

	//sample the data to create the bootstrap samples
	void createBootstraps(int num);
	//clear the bootstraps
	void clearBootstraps();

	//load match score data
	void readScoreData(const std::string filename);

	//calculate the normalized score distribution histograms
	void calculateScoreDistributions(double minScore, double maxScore, int numBins);
	//print the score distributions to an R file
	void printDistsToR(std::string rFilename, std::string label);
	
};

#endif
