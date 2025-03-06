/* coding: utf-8
__version__ = "2.0"
__author__ = "Cristobal Pais"
*/

// Include classes
#include "FBPDataGenerator.h"
#include "FBP5.0.h"
#include "ReadCSV.h"
#include "WriteCSV.h"
#include "ReadArgs.h"

// Include libraries
#include <stdexcept>
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
#include <chrono>

using namespace std;

// Global Variables (DFs with cells and weather info)
inputs * df_ptr;
weatherDF * wdf_ptr;
weatherDF wdf[150];
inputs df [7000000];

/******************************************************************************
															Constructor 
*******************************************************************************/
Cell2Fire::Cell2Fire(arguments _args) : CSVWeather(_args.InFolder + "Weather.csv", ","), 
                                        CSVForest(_args.InFolder + "Forest.asc", " ")
    {
    // Populate arguments from command line into the Cell2Fire object
	this->args = _args;
	this->args_ptr = &this->args;
	
	/********************************************************************
                    Initialize fuel coefficients for FBP
	********************************************************************/
	// fuel_coefs coefs[18];	
	this->coef_ptr = &this->coefs[0];
	setup_const(this->coef_ptr);
	
	/********************************************************************
                        Read Instance from csv files...
	********************************************************************/
	/* Forest DataFrame */
	std::string filename = this->args.InFolder + args.InData;
	std::string sep = ",";
	CSVReader CSVParser(filename, sep);
	
	// Populate DF 
	std::vector<std::vector<std::string>> DF = CSVParser.getData();
	std::cout << "DataFrame from instance " << filename << std::endl;
    if(args.verbose){
        CSVParser.printData(DF);
    }
	
	// Create empty df with 
	df_ptr = & df[0];
                                            
	// Populate the df nrows objects
    int nrows = args.TotalSims;
	CSVParser.parseDF(df_ptr, DF, nrows);
                                            
	// Headers
    std::cout << "" << std::endl;
    std::cout << "HROS" << ',' << "FROS" << ',' << "BROS" << std::endl;                                            
    
    // Calculate parameters
    int i = 0;
    double hros = 0.;
    double fros = 0.;
    double bros = 0.;
                                            
    // Clean structures
    main_outs mainstruct;
    snd_outs sndstruct;
    fire_struc headstruct, backstruct, flankstruct;
    
    // Calculation loop
    for (i = 0; i < nrows; i++) {
        // Pointer to row
        df_ptr = & df[i];
        
        // Weather
        df_ptr->waz = 0.;
        //df_ptr->ws = 10.;
        //df_ptr->ffmc = 90.;
        //df_ptr->bui = 99.;	
        
        // Calculate ROS
        calculate(df_ptr, this->coef_ptr, &mainstruct, &sndstruct, &headstruct, &flankstruct, &backstruct);
       
        // Head, Flank, and Back ROS
        hros = headstruct.ros;
        fros = flankstruct.ros;
        bros = backstruct.ros;
        
        // Print
        std::cout << hros << ',' << fros << ',' << bros << std::endl;
    }
    
    /*  ROSs DEBUG  */
	if(args.verbose){
		std::cout << "*********** ROSs debug ************" << std::endl;
		std::cout <<  "-------Input Structure--------" << std::endl;
		std::cout <<  "fueltype: " << df_ptr->fueltype << std::endl;
		std::cout <<  "ffmc: " << df_ptr->ffmc << std::endl;
		std::cout <<  "ws: " << df_ptr->ws << std::endl;
		std::cout <<  "gfl: " << df_ptr->gfl << std::endl;
		std::cout <<  "bui: " << df_ptr->bui << std::endl;
		std::cout <<  "lat: " << df_ptr->lat << std::endl;
		std::cout <<  "lon: " << df_ptr->lon << std::endl;
		std::cout <<  "time: " << df_ptr->time << std::endl;
		std::cout <<  "pattern: " << df_ptr->pattern << std::endl;
		std::cout <<  "mon: " << df_ptr->mon << std::endl;
		std::cout <<  "jd: " << df_ptr->jd << std::endl;
		std::cout <<  "jd_min: " << df_ptr->jd_min << std::endl;
		std::cout <<  "waz: " << df_ptr->waz << std::endl;
		std::cout <<  "ps: " << df_ptr->ps << std::endl;
		std::cout <<  "saz: " << df_ptr->saz << std::endl;
		std::cout <<  "pc: " << df_ptr->pc << std::endl;
		std::cout <<  "pdf: " << df_ptr->pdf << std::endl;
		std::cout <<  "cur: " << df_ptr->cur << std::endl;
		std::cout <<  "elev: " << df_ptr->elev << std::endl;
		std::cout <<  "hour: " << df_ptr->hour << std::endl;
		std::cout <<  "hourly: " << df_ptr->hourly << std::endl;
		std::cout <<  "\n-------Mainout Structure--------" << std::endl;
		std::cout << "hffmc: " << mainstruct.hffmc << std::endl;
		std::cout << "sfc: " << mainstruct.sfc << std::endl;
		std::cout << "csi: " << mainstruct.csi << std::endl;
		std::cout << "rso: " << mainstruct.rso << std::endl;
		std::cout << "fmc: " << mainstruct.fmc << std::endl;
		std::cout << "sfi: " << mainstruct.sfi << std::endl;
		std::cout << "rss: " << mainstruct.rss << std::endl;
		std::cout << "isi:" << mainstruct.isi << std::endl;
		std::cout << "be:" << mainstruct.be << std::endl;
		std::cout << "sf:" << mainstruct.sf << std::endl;
		std::cout << "raz: " << mainstruct.raz << std::endl;
		std::cout << "wsv:" << mainstruct.wsv << std::endl;
		std::cout << "ff: " << mainstruct.ff << std::endl;
		std::cout << "jd_min:" << mainstruct.jd_min << std::endl;
		std::cout << "jd:" << mainstruct.jd << std::endl;
		std::cout << "covertype: " << mainstruct.covertype << std::endl;
		std::cout <<  "\n-------Headout Structure--------" << std::endl;
		std::cout <<  "ros: " << headstruct.ros << std::endl;
		std::cout <<  "dist: " << headstruct.dist << std::endl;
		std::cout <<  "rost: " << headstruct.rost << std::endl;
		std::cout <<  "cfb: " << headstruct.cfb << std::endl;
		std::cout <<  "fc: " << headstruct.fc << std::endl;
		std::cout <<  "cfc: "<< headstruct.cfc << std::endl;
		std::cout <<  "time: " << headstruct.time << std::endl;
		std::cout <<  "rss: " << headstruct.rss << std::endl;
		std::cout <<  "isi: " << headstruct.isi << std::endl;
		std::cout <<  "fd: " << headstruct.fd << std::endl;
		std::cout <<  "fi: " << headstruct.fi << std::endl;
		std::cout <<  "\n------- Flank Structure--------" << std::endl;
		std::cout <<  "ros: " << flankstruct.ros << std::endl;
		std::cout <<  "dist: " << flankstruct.dist << std::endl;
		std::cout <<  "rost: " << flankstruct.rost << std::endl;
		std::cout <<  "cfb: " << flankstruct.cfb << std::endl;
		std::cout <<  "fc: " << flankstruct.fc << std::endl;
		std::cout <<  "cfc: "<< flankstruct.cfc << std::endl;
		std::cout <<  "time: " << flankstruct.time << std::endl;
		std::cout <<  "rss: " << flankstruct.rss << std::endl;
		std::cout <<  "isi: " << flankstruct.isi << std::endl;
		std::cout <<  "fd: " << flankstruct.fd << std::endl;
		std::cout <<  "fi: " << flankstruct.fi << std::endl;
		std::cout <<  "\n------- Back Structure--------" << std::endl;
		std::cout <<  "ros: " << backstruct.ros << std::endl;
		std::cout <<  "dist: " << backstruct.dist << std::endl;
		std::cout <<  "rost: " << backstruct.rost << std::endl;
		std::cout <<  "cfb: " << backstruct.cfb << std::endl;
		std::cout <<  "fc: " << backstruct.fc << std::endl;
		std::cout <<  "cfc: "<< backstruct.cfc << std::endl;
		std::cout <<  "time: " << backstruct.time << std::endl;
		std::cout <<  "rss: " << backstruct.rss << std::endl;
		std::cout <<  "isi: " << backstruct.isi << std::endl;
		std::cout <<  "fd: " << backstruct.fd << std::endl;
		std::cout <<  "fi: " << backstruct.fi << std::endl;
		std::cout << "*********** ROSs debug ************" << std::endl;
	}
	
}
/*******************************************************************************
                                   Main Program	
*******************************************************************************/
int main(int argc, char * argv[]){
	// Read Arguments
	std::cout << "------ Command line values ------\n";
	arguments args;
	arguments * args_ptr = &args;
	parseArgs(argc, argv, args_ptr);
	
	// Initialize Instance
	Cell2Fire Forest(args);
	
	return 0;
}
