#include "ROC.h"
#include <algorithm>
#include <cmath>
#include <limits>

//default constructor
ROC::ROC(){
	//does nothing
}

ROC::ROC(ScoreData& scores, Parameters& p){
	calculateROC_bootstrap(scores, p);
}

//calculate a single ROC point by point
void ROC::calculateSingleROC(const std::vector<double>& matchScores, const std::vector<double>& nonmatchScores){
	//ensure that the combined score count is acceptable
	std::vector<ROCpoint> tmp;
	if(matchScores.size() + nonmatchScores.size() > tmp.max_size()){
		std::cout<<"Too many scores. Cannot process (VECTOR LIMIT)"<<std::endl;
		exit(-1);
	} 

	//declare the ROC to hold the points
	std::vector<ROCpoint> roc;

	//must sort each of the vectors
	std::vector<double> sortedMatch = matchScores;
	std::vector<double> sortedNonmatch = nonmatchScores;
//	std::cout<<"\t\tSorting the scores"<<std::endl;
	sort(sortedMatch.begin(), sortedMatch.end());
	sort(sortedNonmatch.begin(), sortedNonmatch.end());
	
	//loop through the scores incrementing the TAR and FAR as we go along
	unsigned long long int TA_count = 0; unsigned long long int FA_count = 0;
	std::vector<double>::iterator match_it = sortedMatch.begin();
	std::vector<double>::iterator nonmatch_it = sortedNonmatch.begin();
	ROCpoint tmpPt;
	double prev_val = -1.0;
	double current_val, match_val, nonmatch_val;
	bool TA_flag;
	//flags to know if scores exists
	bool match_flag = true;
	bool nonmatch_flag = true;
	bool roc_flag = true;
	double TAR, FAR;
//	std::cout<<"\t\tCalculating the ROC"<<std::endl;
	while(roc_flag){
		//get a value from match scores if possible
		if(match_flag && match_it != sortedMatch.end()){
			match_val = *match_it;
		}
		else{
			match_flag = false;
			match_val = std::numeric_limits<double>::infinity();
		}

		//get a value from nonmatch scores if possible
		if(nonmatch_flag && nonmatch_it != sortedNonmatch.end()){
			nonmatch_val = *nonmatch_it;
		}
		else{
			nonmatch_flag = false;
			nonmatch_val = std::numeric_limits<double>::infinity();
		}
	
		//break if out of scores
		if((!match_flag) && (!nonmatch_flag)){
			roc_flag = false;
			break;
		}

		//get the next score from match or nonmatch
		if(match_val <= nonmatch_val){
			current_val = match_val;
			++match_it;
			TA_flag = true;
		}
		else{
			current_val = nonmatch_val;
			++nonmatch_it;
			TA_flag = false;
		}

		//check if the current score is the same as the previous score
		if( dbleq(current_val,prev_val) ){
			//do not output an ROC point here
			//but still need to incrememnt counts
			if(TA_flag){
				TA_count++;
			}
			else
				FA_count++;
			
			//just continue to next score
			continue;
		}
		else{
			//need to output an ROC point for the previous values
			TAR = ((double) TA_count) / matchScores.size();
			FAR = ((double) FA_count) / nonmatchScores.size();
			tmpPt.TA = TAR;
			tmpPt.FA = FAR;
			roc.push_back(tmpPt);
			
			//increment the counts
			if(TA_flag)
				TA_count++;
			else
				FA_count++;

			//update prev_val
			prev_val = current_val;
		}
	}

	//output the last point
	TAR = ((double) TA_count) / matchScores.size();
	FAR = ((double) FA_count) / nonmatchScores.size();
	tmpPt.TA = TAR;
	tmpPt.FA = FAR;
	roc.push_back(tmpPt);
	
//	std::cout<<"\t\tPruning the ROC curve"<<std::endl;

	//prune the ROC down to the desired number of points
	double delta_FAR;
	std::vector<ROCpoint>::iterator roc_it;
	roc_it = roc.begin();

	std::vector<ROCpoint> pruned_roc;
	ROCpoint point;
	point.TA = 0.0; point.FA = 0.0;
	pruned_roc.push_back(point);

	double ratio;
	bool start_flag = false;
	
	FAR = 0.000001;
	int ii = 0;
	do{
	
		//increase resolution on left side of curve
		if(FAR < 0.001){
			delta_FAR = 0.000001;
		}
		else if(FAR < 0.01){
			delta_FAR = 0.0001;
		}
		else {
			delta_FAR = 0.001;
		}
		//find the point in the ROC curve with the corresponding FAR
		/*FAR = delta_FAR * ii;
		while(roc_it->FA <= FAR){++roc_it;}
		ratio = (FAR - (roc_it-1)->FA) / (roc_it->FA - (roc_it-1)->FA);
		TAR = (roc_it-1)->TA + ratio * (roc_it->TA - (roc_it-1)->TA);
		pruned_roc[ii].TA = TAR;
		pruned_roc[ii].FA = FAR;
		*/

		while(roc_it->FA < FAR){++roc_it; start_flag = true;}
		if(start_flag){
			ratio = (FAR - (roc_it-1)->FA) / (roc_it->FA - (roc_it-1)->FA);
			TAR = (roc_it-1)->TA + ratio * (roc_it->TA - (roc_it-1)->TA);
		}
		else{
			ratio = (FAR) / (roc_it->FA);
			TAR = ratio * (roc_it->TA);
		}

		point.TA = TAR;
		point.FA = FAR;
		pruned_roc.push_back(point);
		ii++;
		FAR += delta_FAR;
	}while(FAR < 1.0);

	point.TA = 1.0;
	point.FA = 1.0;
	pruned_roc.push_back(point);

	//push the pruned roc on the list
	points.push_back(pruned_roc);
	//points.push_back(roc);
}

//calculate ROCs and bootstrap here
void ROC::calculateROC_bootstrap(ScoreData& scores, Parameters& p){
	//calcualte ROC for the original data
	//calculateSingleROC(scores.match[0],scores.nonmatch[0]);
	if(p.calcOrig){
		std::cout<<"\tOriginal Distribution"<<std::endl;
		calculateSingleROC(scores.match[0],scores.nonmatch[0]);
	}
	
	time_t time1,diff;
	//for each bootstrap
	for(int ii = 1; ii <= p.numBootstraps; ii++){
//		std::cout<<"\tOn bootstrap "<<ii<<std::endl;
		time1 = time(NULL);
		//calculate a bootstrap in the ScoreData object
		//this should clear the old boostraps each time
		scores.createBootstraps(1);
		//calculate a single ROC for this bootstrap
		calculateSingleROC(scores.match[1], scores.nonmatch[1]);
		diff = time(NULL) - time1;
//		std::cout<<"\t\ttook "<<diff<<" seconds"<<std::endl;
	}
}

//calculate the stats for the group of generated ROC curves
void ROC::calculateStats(){
	//needed variables
	double TAmean = 0.0; double TAstddev = 0.0;
	double FAmean = 0.0; double FAstddev = 0.0;
	double TAsum = 0.0; double FAsum = 0.0;
	double TA_CI_95_low = 0.0; double TA_CI_95_high = 0.0;
	double FA_CI_95_low = 0.0; double FA_CI_95_high = 0.0;

	//allocate the stats vector
	stats = std::vector<ROCstat>(points[0].size());

	//allocate a vector to hold the sorted TAs and FAs
	std::vector<double> sortedTAs, sortedFAs;

	//for each threshold, calculate the mean and stddev of the TA and FA
	for(int ii = 0; ii < points[0].size(); ii++){
		//clear the sorted vectors
		sortedTAs.erase(sortedTAs.begin(),sortedTAs.end());
		sortedFAs.erase(sortedFAs.begin(),sortedFAs.end());

		//calculate the mean
		TAsum = 0.0; FAsum = 0.0;
		for(int jj = 0; jj < points.size(); jj++){
			TAsum += points[jj][ii].TA;
			FAsum += points[jj][ii].FA;
			sortedTAs.push_back(points[jj][ii].TA);
			sortedFAs.push_back(points[jj][ii].FA);
		}
		TAmean = TAsum / points.size();
		FAmean = FAsum / points.size();

		//sort the TA and FA vectors
		sort(sortedTAs.begin(),sortedTAs.end());
		sort(sortedFAs.begin(),sortedFAs.end());
		//calculate the 2.5% and 97.5% interval based on histogram
		int indices[2];
		indices[0] = (int) (floor( 0.025 * (points.size()-1)) + 0.5); 	
		indices[1] = (int) (ceil( 0.975 * (points.size()-1)) + 0.5); 	
		TA_CI_95_low = sortedTAs[indices[0]];
		TA_CI_95_high = sortedTAs[indices[1]];
		FA_CI_95_low = sortedFAs[indices[0]];
		FA_CI_95_high = sortedFAs[indices[1]];
		//write to the temp object and then push on vector
		CI temp;
		temp.TA_low = TA_CI_95_low;
		temp.TA_high = TA_CI_95_high;
		temp.FA_low = FA_CI_95_low;
		temp.FA_high = FA_CI_95_high;
		CIs.push_back(temp);

		//calculate the standard deviation
		TAsum = 0.0; FAsum = 0.0;
		TAstddev = 0.0; FAstddev = 0.0;
		for(int jj = 0; jj < points.size(); jj++){
			TAsum += (points[jj][ii].TA - TAmean) * (points[jj][ii].TA - TAmean);
			FAsum += (points[jj][ii].FA - FAmean) * (points[jj][ii].FA - FAmean);
		}
		TAstddev = sqrt(TAsum/(points.size()-1)); 
		FAstddev = sqrt(FAsum/(points.size()-1)); 

		//put the values in the stat vector
		stats[ii].TAmean = TAmean;
		stats[ii].FAmean = FAmean;
		stats[ii].TAstd = TAstddev;
		stats[ii].FAstd = FAstddev;
	}
}

//print the statistics
void ROC::printStats(){
	for(int ii = 0; ii < stats.size(); ii++){
		std::cout<<"Threshold : "<<thresholds[ii]<<std::endl;
		std::cout<<"\tMeans (TA,FA) : "<<stats[ii].TAmean<<"\t"<<stats[ii].FAmean<<std::endl;
		std::cout<<"\tStdDev(TA,FA) : "<<stats[ii].TAstd<<"\t"<<stats[ii].FAstd<<std::endl;
	}
}

//print in R format to a file
void ROC::printStatsToR(Parameters& p){
	//open the file
	std::ofstream rFile(p.rocOut.c_str(), std::ofstream::out); 	
	if(!rFile.good()){
		std::cout<<"Unable to open file"<<std::endl;
		return;
	}
	
	bool first = true;

	//write the TA means
	rFile<<"meanTA_"<<p.label<<"=[";
	for(int ii = 0; ii < stats.size(); ii++){
		if(first) {
			rFile<<stats[ii].TAmean;
			first =  false;
		}
		else{
			rFile<<","<<stats[ii].TAmean;
		}
	}	
	rFile<<"];"<<std::endl;
	first = true;
	
	//write the FA means
	rFile<<"meanFA_"<<p.label<<"=[";
	for(int ii = 0; ii < stats.size(); ii++){
		if(first) {
			rFile<<stats[ii].FAmean;
			first =  false;
		}
		else{
			rFile<<","<<stats[ii].FAmean;
		}
	}	
	rFile<<"];"<<std::endl;
	first = true;

	//if bootstrapped then calculate the stddevs
	if(points.size() > 1){	
		//write the TA stddev
		rFile<<"stdTA_"<<p.label<<"=[";
		for(int ii = 0; ii < stats.size(); ii++){
			if(first) {
				rFile<<stats[ii].TAstd;
				first =  false;
			}
			else{
				rFile<<","<<stats[ii].TAstd;
			}
		}	
		rFile<<"];"<<std::endl;
		first = true;

		//write the FA stddev
		rFile<<"stdFA_"<<p.label<<"=[";
		for(int ii = 0; ii < stats.size(); ii++){
			if(first) {
				rFile<<stats[ii].FAstd;
				first =  false;
			}
			else{
				rFile<<","<<stats[ii].FAstd;
			}
		}	
		rFile<<"];"<<std::endl;
		first = true;

		//write the TA CI_low
		rFile<<"CITA_low_"<<p.label<<"=[";
		for(int ii = 0; ii < CIs.size(); ii++){
			if(first) {
				rFile<<CIs[ii].TA_low;
				first =  false;
			}
			else{
				rFile<<","<<CIs[ii].TA_low;
			}
		}	
		rFile<<"];"<<std::endl;
		first = true;

		//write the TA CI_high
		rFile<<"CITA_high_"<<p.label<<"=[";
		for(int ii = 0; ii < CIs.size(); ii++){
			if(first) {
				rFile<<CIs[ii].TA_high;
				first =  false;
			}
			else{
				rFile<<","<<CIs[ii].TA_high;
			}
		}	
		rFile<<"];"<<std::endl;
		first = true;

		//write the FA CI_low
		rFile<<"CIFA_low_"<<p.label<<"=[";
		for(int ii = 0; ii < CIs.size(); ii++){
			if(first) {
				rFile<<CIs[ii].FA_low;
				first =  false;
			}
			else{
				rFile<<","<<CIs[ii].FA_low;
			}
		}	
		rFile<<"];"<<std::endl;
		first = true;

		//write the FA CI_high
		rFile<<"CIFA_high_"<<p.label<<"=[";
		for(int ii = 0; ii < CIs.size(); ii++){
			if(first) {
				rFile<<CIs[ii].FA_high;
				first =  false;
			}
			else{
				rFile<<","<<CIs[ii].FA_high;
			}
		}	
		rFile<<"];"<<std::endl;
		first = true;
	}

	//print plotting information to plot the roc
//	rFile<<"pdf(\"roc_"<<p.label<<".pdf\");"<<std::endl;
//	rFile<<"plot(meanFA_"<<p.label<<",meanTA_"<<p.label<<",type='l');"<<std::endl;
//	rFile<<"dev.off();"<<std::endl;
	
	rFile.close();

}


//function to compare doubles
//returns true if they are equal, false otherwise
inline bool ROC::dbleq(double a, double b){
	if( ((b - std::numeric_limits<double>::epsilon()) < a) && (a < (b + std::numeric_limits<double>::epsilon())) ){
		return true;
	}
	return false;
}
