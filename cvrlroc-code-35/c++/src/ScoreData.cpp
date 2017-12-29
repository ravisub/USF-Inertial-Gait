#include "ScoreData.h"
#include "MersenneTwister.h"
#include <limits>
#include <algorithm>

ScoreData::ScoreData(){
	numBoot = 0;
}

ScoreData::ScoreData(const std::string filename){
	readScoreData(filename);
	RNG = MTRand(1);
}

//copy constructor
ScoreData::ScoreData(const ScoreData& sd){
	match = sd.match;
	nonmatch = sd.nonmatch;
	matchDist = sd.matchDist;
	nonmatchDist = sd.nonmatchDist;
	binMax = sd.binMax;
}

//convert higher is better to lower is better
void ScoreData::convertHigherToLower(){
	//find the max element in the base scores
	double maxMatch = *std::max_element(match[0].begin(), match[0].end());
	double maxNonmatch = *std::max_element(nonmatch[0].begin(), nonmatch[0].end());
	double max;
	if(maxMatch >= maxNonmatch)
		max = maxMatch;
	else
		max = maxNonmatch;

	//loop through and convert every score to max - score
	std::vector<double>::iterator it;
	for(it = match[0].begin(); it != match[0].end(); ++it){
		*it = max - *it;
	}
	for(it = nonmatch[0].begin(); it != nonmatch[0].end(); ++it){
		*it = max - *it;
	}
}

//create bootstrap samples of the data
void ScoreData::createBootstraps(int num){
	//if the score vectors already contain bootstraps, delete them
	clearBootstraps();

	//push num vectors on the score vectors
	for(int ii = 0; ii < num; ii++){
		match.push_back(std::vector<double>(match[0].size()));
		nonmatch.push_back(std::vector<double>(nonmatch[0].size()));
	}

	//check if the number of match or nonmatch scores is greater than limit of unsigned long int
	if(match[0].size() > std::numeric_limits<unsigned long int>::max() | 
		nonmatch[0].size() > std::numeric_limits<unsigned long int>::max()){
		std::cout<<"Distributions larger than "<<std::numeric_limits<unsigned long int>::max()<<" are not currently supported"<<std::endl;
		exit(-1);
	}

	//if we made it this far, the distribution is smaller than the unsigned long int limit
	//set up the RNG
//	MTRand RNG(time(NULL));
//	MTRand RNG(1);
	
	//loop through each of the distributions and bootstrap
	for(int ii = 1; ii <= num; ii++){
		std::vector<double>::iterator it;
		for(it = match[ii].begin(); it != match[ii].end(); ++it){
			//generate a random number and place it at the current position
			*it = match[0][RNG.randInt(match[0].size())];
		}
		//do the same for nonmatch
		for(it = nonmatch[ii].begin(); it != nonmatch[ii].end(); ++it){
			//generate a random number and place it at the current position
			*it = nonmatch[0][RNG.randInt(nonmatch[0].size())];
		}
	}
}

//clear the bootstraps
void ScoreData::clearBootstraps(){
	if(match.size() > 1){
		//std::cout<<"Clearing match bootstraps"<<std::endl;
		match.erase(match.begin()+1,match.end());
	}
	if(nonmatch.size() > 1){
		//std::cout<<"Clearing nonmatch bootstraps"<<std::endl;
		nonmatch.erase(nonmatch.begin()+1,nonmatch.end());
	}
	numBoot = 0;
}

void ScoreData::readScoreData(const std::string filename){
	//open the file
	std::ifstream scoreFile(filename.c_str(), std::ifstream::in);

	//clear the vectors
	match.erase(match.begin(), match.end());
	nonmatch.erase(nonmatch.begin(), nonmatch.end());
	//push a new row onto each vector
	match.push_back(std::vector<double>());
	nonmatch.push_back(std::vector<double>());
	
	//go line by line and extract the scores
	double score = -1.0; int result = -1;
	while(scoreFile.good()){
		scoreFile>>result>>score;
		if(result == 1){
			//match
			match[0].push_back(score); 
		}
		else if (result == 0){
			//nonmatch
			nonmatch[0].push_back(score);
		}
	}
}

//calculate the normalized score distribution histograms
void ScoreData::calculateScoreDistributions(double minScore, double maxScore, int numBins){
	//declare a vector to store the binMaximums
	double delta = (maxScore - minScore) / numBins;
	for(double score = minScore; score < maxScore - delta; score += delta){
		binMax.push_back(score + delta);
	}

	matchDist = std::vector<double>(binMax.size(),0.0);
	nonmatchDist = std::vector<double>(binMax.size(),0.0);
	//for each score, find the correct bin
	std::vector<double>::iterator score_iter, bin_size_iter, bin_iter;
	for(score_iter = match[0].begin(); score_iter != match[0].end(); ++score_iter){
		//find the first bin it fits into
		bin_iter = matchDist.begin();
		for(bin_size_iter = binMax.begin(); bin_size_iter != binMax.end(); ++bin_size_iter){
			if( *score_iter <= *bin_size_iter){
				(*bin_iter) += 1.0;
				break;
			}
			++bin_iter;
		}
	}
	for(score_iter = nonmatch[0].begin(); score_iter != nonmatch[0].end(); ++score_iter){
		bin_iter = nonmatchDist.begin();
		//find the first bin it fits into
		for(bin_size_iter = binMax.begin(); bin_size_iter != binMax.end(); ++bin_size_iter){
			if( *score_iter <= *bin_size_iter){
				(*bin_iter) += 1.0;
				break;
			}
			++bin_iter;
		}
	}

	//normalize the bin sizes
	for(bin_iter = matchDist.begin(); bin_iter != matchDist.end(); ++bin_iter){
		*bin_iter = (*bin_iter) / match[0].size();
	}
	for(bin_iter = nonmatchDist.begin(); bin_iter != nonmatchDist.end(); ++bin_iter){
		*bin_iter = (*bin_iter) / nonmatch[0].size();
	}
}

//print the score distributions to an R file
void ScoreData::printDistsToR(std::string rFilename, std::string label){
	//open the file
	std::ofstream rFile(rFilename.c_str(), std::ofstream::out); 	
	if(!rFile.good()){
		std::cout<<"Unable to open file"<<std::endl;
		return;
	}
	
	bool first = true;
	std::vector<double>::iterator it;
	
	//write the bin centers
	double delta = binMax[1] - binMax[0];
	double center;
	rFile<<"binCenters_"<<label<<"<-c(";	
	for(it = binMax.begin(); it != binMax.end(); ++it){
		if(first) {
			center = *it - delta;
			rFile<<center;
			first =  false;
		}
		else{
			center = *it - delta;
			rFile<<","<<center;
		}
	}	
	rFile<<")"<<std::endl;
	first = true;

	//write the matchDist
	rFile<<"matchDist_"<<label<<"<-c(";	
	for(it = matchDist.begin(); it != matchDist.end(); ++it){
		if(first) {
			rFile<<*it;
			first =  false;
		}
		else{
			rFile<<","<<*it;
		}
	}	
	rFile<<")"<<std::endl;
	first = true;
	
	//write the nonmatchDist
	rFile<<"nonmatchDist_"<<label<<"<-c(";	
	for(it = nonmatchDist.begin(); it != nonmatchDist.end(); ++it){
		if(first) {
			rFile<<*it;
			first =  false;
		}
		else{
			rFile<<","<<*it;
		}
	}	
	rFile<<")"<<std::endl;
	first = true;

	rFile.close();
}
