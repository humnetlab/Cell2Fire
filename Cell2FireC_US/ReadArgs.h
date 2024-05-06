// Definition
#ifndef READARGS
#define READARGS

// Include libraries
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <algorithm>

/*
*   Args structure
*/
typedef struct{ 
    std::string InFolder, OutFolder, WeatherOpt, HarvestPlan, KOption;
    bool OutMessages, Trajectories, NoOutput, verbose, Ignitions, OutputGrids, FinalGrid, PromTuned, Stats, BBOTuning, AllowCROS, CenteredDistribution, NoAngleOffset;
    float ROSCV, ROSThreshold, CROSThreshold, HFIThreshold, HFactor, FFactor, BFactor, EFactor, FirePeriodLen;
    float CBDFactor, CCFFactor, ROS10Factor, CROSActThreshold;
    float KFactor1, KFactor2, KFactor3, KFactor4, KFactor5;
    int MinutesPerWP, MaxFirePeriods, TotalYears, TotalSims, NWeatherFiles, IgnitionRadius, EllipticalOption, LBFormula, seed;
    std::unordered_set<int>  HCells, BCells;
} arguments;


char* getCmdOption(char ** begin, char ** end, const std::string & option);

bool cmdOptionExists(char** begin, char** end, const std::string& option);

void parseArgs(int argc, char * argv[], arguments * args_ptr);

void printArgs(arguments args);


#endif