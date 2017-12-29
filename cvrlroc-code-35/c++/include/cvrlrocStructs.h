#ifndef CVRLROCSTRUCTS_H
#define CVRLROCSTRUCTS_H

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <ctime>
#include <vector>
#include <fstream>

//represents a point on a ROC curve
struct ROCpoint{
	double TA;
	double FA;
};

//represents the statistics for set of bootstrapped points
struct ROCstat{
	double TAmean;
	double TAstd;
	double FAmean;
	double FAstd;
};

//contains CI information
struct CI{
	double TA_low;
	double TA_high;
	double FA_low;
	double FA_high;
};

#endif
