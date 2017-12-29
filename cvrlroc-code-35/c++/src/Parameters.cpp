#include "Parameters.h"
#include <iostream>
#include <cstdlib>

Parameters::Parameters() : processed(false) {}

Parameters::Parameters(int argc, char** argv){
	processCommandLine(argc, argv);
	processed = true;
}

//display the usage information
void Parameters::usage(){
	std::cout<<"USAGE : ./CVRLROC [options]"<<std::endl
		<<"Required options are noted"<<std::endl<<std::endl
		<<"-s, --scoreFile  [scoreFile]       : name of the input file containing the score (REQUIRED)"<<std::endl
		<<"-r, --rocFile    [rocFile]         : name of the output file to store the ROC (REQUIRED)"<<std::endl<<std::endl
		<<"-n, --numBoot    [numBootstraps]   : number of bootstraps to use (DEFAULT = 0)"<<std::endl
		<<"-l, --label      [label]           : will append _[label] to vector names in output files (DEFAULT = \"\")"<<std::endl
		<<"-b, --highBettter                  : inform the program that higher scores are better (DEFAULT OFF)"<<std::endl
		<<"-d, --calcDists  [distFile]        : if set, calculate the distribution statistics (DEFAULT OFF)"<<std::endl<<std::endl
		<<"-c, --calcOriginal                 : if set, calculates roc for original distribution (DEFAULT OFF if bootstrapping on,"<<std::endl
		<<"                                                                                               ON if no bootstrapping)"<<std::endl
		<<"-h, --help                         : display the usage information"<<std::endl<<std::endl;
	
	exit(-1);
}

//process the command line arguments
void Parameters::processCommandLine(int argc, char** argv){
	//set the default parameters
	numBootstraps = 0;
	label = std::string("");
	calcDists = false;
	distOut = std::string("");
	highBetter = false;
	calcOrig = false;

	//initialize the GetOpt_pp object
	GetOpt::GetOpt_pp ops(argc,argv);

	//check for the help flag
	if(ops >> GetOpt::OptionPresent('h',"help"))
		usage();

	//check for the scoreFile flag
	if(ops >> GetOpt::Option('s',"scoreFile",scoreFile)){}
	else usage();

	//check for the ROC output file
	if(ops >> GetOpt::Option('r',"rocFile",rocOut)){}
	else usage();

	//check for the calcDists flag
	if(ops >> GetOpt::Option('d',"calcDists",distOut))
		calcDists = true;

	//check for the numBoot flag
	if(ops >> GetOpt::Option('n',"numBoot",numBootstraps)){}

	//check for the label flag
	if(ops >> GetOpt::Option('l',"label",label)){}

	//check for the highBetter flag
	if(ops >> GetOpt::OptionPresent('b', "highBetter"))
		highBetter = true;

	//check for the calc orig flag
	if(ops >> GetOpt::OptionPresent('c', "calcOriginal"))
		calcOrig = true;
	
	//set calcOriginal to default if numBoot == 0 and no flag present
	if(!calcOrig && numBootstraps == 0)
		calcOrig = true;
}
