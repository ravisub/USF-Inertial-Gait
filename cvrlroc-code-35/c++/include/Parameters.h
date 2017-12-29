#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "getopt_pp.h"
#include <string>

class Parameters{
public:
	//default constructor
	Parameters();
	//initialize from command line
	Parameters(int argc, char** argv);

	//display usage information
	void usage();

	//process the command line arguments
	void processCommandLine(int argc, char** argv);
	
	//variables for parameter values
	std::string scoreFile;
	std::string rocOut;
	std::string distOut;
	std::string label;
	int numBootstraps;
	bool calcDists;
	bool processed;
	bool highBetter;
	bool calcOrig;
	
};

#endif
