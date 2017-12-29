#include "Parameters.h"
#include "ScoreData.h"
#include "ROC.h"
#include <iostream>

using namespace std;

int main(int argc, char** argv){
	Parameters p(argc,argv);

	//load the scores
	cout<<"Loading the score file"<<endl;
	ScoreData scores(p.scoreFile);

	//if higher scores are better, convert
	if(p.highBetter){
		cout<<"Converting scores to Lower Is Better format"<<endl;
		scores.convertHigherToLower();
	}

	//calculate the ROCs
	cout<<"Calculating ROCs"<<endl;
	ROC roc(scores,p);
	roc.calculateStats();

	//calculate the distribution information if desired
	if(p.calcDists){

	}

	//print the ROC output
	roc.printStatsToR(p);
}	
