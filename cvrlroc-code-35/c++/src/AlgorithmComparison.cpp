#include "AlgorithmComparison.h"
#include <string>

//default constructor
AlgorithmComparison::AlgorithmComparison(){
	//don't do anything
}
//constructor that accepts score data and calculates everything
AlgorithmComparison::AlgorithmComparison(ScoreData& alg1, ScoreData& alg2, int numBoot, int numROCPoint, double FApoint){
	//set the target FA
	FA_target_rate = FApoint;

	//set the pointers for the scores
	a1 = &alg1; a2 = &alg2;
	
	//clear any bootstraps the might be present
	a1->clearBootstraps(); a2->clearBootstraps();

	//calculate the ROCs from the initial data
	roc1 = new ROC(*a1, 0.0, 0.5, numROCPoint);	
	roc2 = new ROC(*a2, 0.0, 0.5, numROCPoint);	

	time_t time1, diff;
	//calculate the ROCs for the bootstrap
	for(int ii = 0; ii < numBoot; ii++){
//		std::cout<<"On bootstrap "<<ii<<std::endl;
		time1 = time(NULL);
		//clear any bootstraps the might be present
		a1->clearBootstraps(); a2->clearBootstraps();
		createSynchronizedBootstraps(2);
		//calculate a single ROC for this bootstrap
		roc1->calculateSingleROC(a1->match[1], a1->nonmatch[1]);
		roc2->calculateSingleROC(a2->match[1], a2->nonmatch[1]);
		diff = time(NULL) - time1;
//		std::cout<<"\t\ttook "<<diff<<" seconds"<<std::endl;
		
	}

	//compute the correlation coefficient for TA @ FA = 0.001
	calculateCorrelationTAatFA(FApoint);
	
	//calculate the Z statistic between the two algorithms for TA @ FA = 0.001
	calculateZStatisticTAatFA();
}

//create bootstraps of corresponding scores
bool AlgorithmComparison::createSynchronizedBootstraps(int num){
	
	//before anything, check to make sure dists are the same size
	if(a1->match[0].size() != a2->match[0].size() ||
		a1->nonmatch[0].size() != a2->nonmatch[0].size()){
		//not same size RETURN
		std::cout<<"Synchronized bootstrap not performed. ScoreDatas are different sizes"<<std::endl;
		return false;
	}

	//clear the bootstraps
	a1->clearBootstraps();
	a2->clearBootstraps();

	//need to account for score sizes bigger than the random number generator 
	//seed the number generator	
//	srand(time(NULL));
//	srand(1);

	//we need to know the max score size in order to know how to generate the random number
	GeneratorIntervals matchIntervals, nonmatchIntervals;
	matchIntervals.numIntervals = 1;
	nonmatchIntervals.numIntervals = 1;
	unsigned long long int match_count = (unsigned long long int) a1->match[0].size();	
	unsigned long long int nonmatch_count = (unsigned long long int) a1->nonmatch[0].size();	
	
	//determine how many intervals must be used for the random number generator
	while(match_count > RAND_MAX){
		matchIntervals.numIntervals++;
		match_count -= RAND_MAX;
	}
	while(nonmatch_count > RAND_MAX){
		nonmatchIntervals.numIntervals++;
		nonmatch_count -= RAND_MAX;
	}
	match_count = (unsigned long long int) a1->match[0].size();	
	nonmatch_count = (unsigned long long int) a2->nonmatch[0].size();	

	//calculate the interval size and the overflow on the last interval
	matchIntervals.intervalSize = match_count / matchIntervals.numIntervals;
	matchIntervals.overflow = match_count - (matchIntervals.intervalSize * matchIntervals.numIntervals);
	nonmatchIntervals.intervalSize = nonmatch_count / nonmatchIntervals.numIntervals;
	nonmatchIntervals.overflow = nonmatch_count - (nonmatchIntervals.intervalSize * nonmatchIntervals.numIntervals);

	//push num - 1 vectors on the score vectors
	for(int ii = 1; ii < num; ii++){
		a1->match.push_back(std::vector<double>(a1->match[0].size()));
		a1->nonmatch.push_back(std::vector<double>(a1->nonmatch[0].size()));
		a2->match.push_back(std::vector<double>(a2->match[0].size()));
		a2->nonmatch.push_back(std::vector<double>(a2->nonmatch[0].size()));
	}

	int interval;
	unsigned long long int position;
	for(int ii = 1; ii < num; ii++){
		//loop through the score vectors and shuffle the elements
		std::vector<double>::iterator a1_score_it, a2_score_it;
		a2_score_it = a2->match[ii].begin();
		for(a1_score_it = a1->match[ii].begin(); a1_score_it != a1->match[ii].end(); ++a1_score_it){
			//generate the random number for the matches
			interval = rand() % matchIntervals.numIntervals;
			if(interval = matchIntervals.numIntervals-1){
				//in the last interval, need to account for overflow
				position = rand() % (matchIntervals.intervalSize + matchIntervals.overflow);
			}
			else {
				position = rand() % matchIntervals.intervalSize;
			}

			*a1_score_it = a1->match[0][interval * matchIntervals.intervalSize + position];
			*a2_score_it = a2->match[0][interval * matchIntervals.intervalSize + position];
			
			++a2_score_it;
		}		
		a2_score_it = a2->nonmatch[ii].begin();
		for(a1_score_it = a1->nonmatch[ii].begin(); a1_score_it != a1->nonmatch[ii].end(); ++a1_score_it){
			//generate the random number for the nonmatches
			interval = rand() % nonmatchIntervals.numIntervals;
			if(interval = nonmatchIntervals.numIntervals-1){
				//in the last interval, need to account for overflow
				position = rand() % (nonmatchIntervals.intervalSize + nonmatchIntervals.overflow);
			}
			else {
				position = rand() % nonmatchIntervals.intervalSize;
			}

			*a1_score_it = a1->nonmatch[0][interval * nonmatchIntervals.intervalSize + position];
			*a2_score_it = a2->nonmatch[0][interval * nonmatchIntervals.intervalSize + position];
			
			++a2_score_it;
		}		

	}
	a1->numBoot = num;
	a2->numBoot = num;
	return true;
}

//calculate the correlateion for TA @ FA
void AlgorithmComparison::calculateCorrelationTAatFA(double FApoint){
	//declare a vector to hold the TA statistics
	std::vector<double> a1_TAs;
	std::vector<double> a2_TAs;
	
	//loop through each ROC and find the point @ FA
	std::vector<std::vector<ROCpoint> >::iterator ROC_iter;
	std::vector<ROCpoint>::iterator point_iter;
	double interpolateFactor, TApoint;
	//do it for roc1
	for(ROC_iter = roc1->points.begin(); ROC_iter != roc1->points.end(); ++ROC_iter){
		for(point_iter = ROC_iter->begin(); point_iter != ROC_iter->end(); ++point_iter){
			if(point_iter->FA < FApoint){
				continue;
			}
			//we have found the first point with FA greater than the given threshold
			//interpolate to get the TA value at the FApoint
			interpolateFactor = (FApoint - (point_iter-1)->FA) / (point_iter->FA - (point_iter-1)->FA);
			TApoint = (point_iter-1)->TA + interpolateFactor * (point_iter->TA - (point_iter-1)->TA);
			a1_TAs.push_back(TApoint);
			break;
		}
	}
	//do it for roc2
	for(ROC_iter = roc2->points.begin(); ROC_iter != roc2->points.end(); ++ROC_iter){
		for(point_iter = ROC_iter->begin(); point_iter != ROC_iter->end(); ++point_iter){
			if(point_iter->FA < FApoint){
				continue;
			}
			//we have found the first point with FA greater than the given threshold
			//interpolate to get the TA value at the FApoint
			interpolateFactor = (FApoint - (point_iter-1)->FA) / (point_iter->FA - (point_iter-1)->FA);
			TApoint = (point_iter-1)->TA + interpolateFactor * (point_iter->TA - (point_iter-1)->TA);
			a2_TAs.push_back(TApoint);
			break;
		}
	}


	//now we need the mean and stddev of the the TAs
	a1_TAmean = 0.0;
	a2_TAmean = 0.0;
	a1_TAstd = 0.0; 
	a2_TAstd = 0.0;

	std::vector<double>::iterator TA_iter;
	//for a1
	for(TA_iter = a1_TAs.begin(); TA_iter != a1_TAs.end(); ++TA_iter){
		a1_TAmean += *TA_iter;
	}
	a1_TAmean = a1_TAmean / a1_TAs.size();
	//now for a2
	for(TA_iter = a2_TAs.begin(); TA_iter != a2_TAs.end(); ++TA_iter){
		a2_TAmean += *TA_iter;
	}
	a2_TAmean = a2_TAmean / a2_TAs.size();

	//compute the standard deviations
	for(TA_iter = a1_TAs.begin(); TA_iter != a1_TAs.end(); ++TA_iter){
		a1_TAstd += (*TA_iter - a1_TAmean) * (*TA_iter - a1_TAmean);
	}
	a1_TAstd = sqrt(a1_TAstd / (a1_TAs.size()-1));
	for(TA_iter = a2_TAs.begin(); TA_iter != a2_TAs.end(); ++TA_iter){
		a2_TAstd += (*TA_iter - a2_TAmean) * (*TA_iter - a2_TAmean);
	}
	a2_TAstd = sqrt(a2_TAstd / (a2_TAs.size()-1));
	
	//now compute the correlation coefficient
	std::vector<double>::iterator a1_iter, a2_iter;
	double a1Sum = 0.0;
	double a1Sum2 = 0.0;
	double a2Sum = 0.0;
	double a2Sum2 = 0.0;
	double a1a2Sum = 0.0;
	int n = a1_TAs.size();

	a2_iter = a2_TAs.begin();
	for(a1_iter = a1_TAs.begin(); a1_iter != a1_TAs.end(); ++a1_iter){
		a1Sum += *a1_iter;
		a1Sum2 += *a1_iter * *a1_iter;
		a2Sum += *a2_iter;
		a2Sum2 += *a2_iter * *a2_iter;
		a1a2Sum += *a1_iter * *a2_iter;
		++a2_iter;
	}
	cc = ( (n * a1a2Sum) - (a1Sum * a2Sum) ) / ( sqrt( (n * a1Sum2) - (a1Sum*a1Sum) ) * sqrt( (n * a2Sum2) - (a2Sum*a2Sum) ) ); 
	
}

//calculate the Z statistic between the two algorithms for TA @ FA
void AlgorithmComparison::calculateZStatisticTAatFA(){
	zStat = (a1_TAmean - a2_TAmean) / sqrt( a1_TAstd*a1_TAstd + a2_TAstd*a2_TAstd - 2*cc*a1_TAstd*a2_TAstd);
}

//print the correlation and z-stat to a file
void AlgorithmComparison::printStatsToR(std::string filename, std::string label1, std::string label2){
	std::ofstream rFile(filename.c_str(), std::ofstream::out); 	
	if(!rFile.good()){
		std::cout<<"Unable to open file"<<std::endl;
		return;
	}
	rFile<<"correlationCoeff_"<<label1<<"_"<<label2<<" <- "<<cc<<std::endl;
	rFile<<"zStat_"<<label1<<"_"<<label2<<" <- "<<zStat<<std::endl;
	rFile<<"FArate_"<<label1<<"_"<<label2<<" <- "<<FA_target_rate<<std::endl;
	rFile<<"TAatFA_"<<label1<<" <- "<<a1_TAmean<<std::endl;
	rFile<<"TAatFA_"<<label2<<" <- "<<a2_TAmean<<std::endl;
	rFile<<"TAatFAstd_"<<label1<<" <- "<<a1_TAstd<<std::endl;
	rFile<<"TAatFAstd_"<<label2<<" <- "<<a2_TAstd<<std::endl;
}
