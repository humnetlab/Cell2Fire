#ifndef FBPDataGenerator
#define FBPDataGenerator 

// Headers
#include "FBP5.0.h"
#include "ReadCSV.h"
#include "WriteCSV.h"
#include "ReadArgs.h"

// Include libraries
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <math.h>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <string.h>
#include <random>
#include <algorithm> 

using namespace std;

class Cell2Fire {
	// DFs 
	private:
		CSVReader CSVWeather;
		CSVReader CSVForest;
			
    public: 
		 // Main inputs
		 arguments args;
		 arguments* args_ptr;
		 fuel_coefs* coef_ptr;
		 fuel_coefs coefs[18];	
		 
		 // Integers
		 int rows = 0;
		 int cols = 0;
		 int weatherPeriod = 0;
		 int year = 1;
		 int gridNumber = 0;
		 int weatherperiod = 0;
			
    
		// Strings	
		 string gridFolder;
	
		 // Vectors
		 std::vector<std::vector<std::string>> DF;
		 std::vector<std::vector<std::string>> WeatherDF;
		 
		 // Constructor
        Cell2Fire(arguments args);
};

#endif
