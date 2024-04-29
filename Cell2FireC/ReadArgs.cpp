// Inclusions
#include "ReadArgs.h"

// Include libraries
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <string>
#include <algorithm>

// Read CMD options
char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

// Check existence
bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}

// Parse options/args
void parseArgs(int argc, char * argv[], arguments * args_ptr)
{
    // Help
    if(cmdOptionExists(argv, argv+argc, "-h")){
        printf("-------------------------------------------\n         Help manual!!! \n-------------------------------------------\n");
    }

    // Empty default
    char empty = '\0';

    // Strings
    //--input-instance-folder
    char * input_folder = getCmdOption(argv, argv + argc, "--input-instance-folder");
    if (input_folder){
        printf("InFolder: %s \n", input_folder);
    }
    else input_folder = &empty;

    //--output-folder
    char * output_folder = getCmdOption(argv, argv + argc, "--output-folder");
    if (output_folder){
        printf("OutFolder: %s \n", output_folder );
    }
    else output_folder = &empty;

    //--weather
    char * input_weather = getCmdOption(argv, argv + argc, "--weather");
    if (input_weather){
        printf("WeatherOpt: %s \n", input_weather);
    }
    else input_weather = &empty;

    //--HarvestPlan
    char * input_hplan = getCmdOption(argv, argv + argc, "--HarvestPlan");
    if (input_hplan){
        printf("HarvestPlan: %s \n", input_hplan);
    }
    else input_hplan = &empty;

    //--kopt (MHK)
    char * kopt = getCmdOption(argv, argv + argc, "--kopt");
    if (kopt){
        printf("Optimal K factors: %s \n", kopt);
    }
    else kopt = &empty;


    // Booleans
    bool out_messages = false;
    bool out_trajectories = false;
    bool no_output = false; 
    bool verbose_input = false;
    bool input_ignitions = false;
    bool out_grids = false;
    bool out_finalgrid = false;
    bool prom_tuned = false;
    bool out_stats = false;	
    bool bbo_tuning = false;
    bool allow_cros = false;
    bool centered_distribution = false;
    bool deactivate_offset = false;
    // bool kopt = false;

    //--out-messages
    if(cmdOptionExists(argv, argv+argc, "--output-messages")){
        out_messages = true;
        printf("OutMessages: %d \n", out_messages);
    }

    //--trajectories
    if(cmdOptionExists(argv, argv+argc, "--trajectories")){
        out_trajectories = true;
        printf("Trajectories: %d \n", out_trajectories);
    }

    //--no-output
    if(cmdOptionExists(argv, argv+argc, "--no-output")){
        no_output = true;
        printf("noOutput: %d \n", no_output);
    }

    //--verbose
    if(cmdOptionExists(argv, argv+argc, "--verbose")){
        verbose_input = true;
        printf("verbose: %d \n", verbose_input);
    }

    //--ignitions
    if(cmdOptionExists(argv, argv+argc, "--ignitions")){
        input_ignitions = true;
        printf("Ignitions: %d \n", input_ignitions);
    }

    //--grids
    if(cmdOptionExists(argv, argv+argc, "--grids")){
        out_grids = true;
        printf("OutputGrids: %d \n", out_grids);
    }

    //--final-grid
    if(cmdOptionExists(argv, argv+argc, "--final-grid")){
        out_finalgrid = true;
        printf("FinalGrid: %d \n", out_finalgrid);
    }

    //--Prom_tuned
    if(cmdOptionExists(argv, argv+argc, "--PromTuned")){
        prom_tuned = true;
        printf("PromTuned: %d \n", prom_tuned);
    }

    //--statistics
    if(cmdOptionExists(argv, argv+argc, "--statistics")){
        out_stats = true;
        printf("Statistics: %d \n", out_stats);
    }

    //--bbo
    if(cmdOptionExists(argv, argv+argc, "--bbo")){
        bbo_tuning = true;
        printf("BBOTuning: %d \n", out_stats);
    }

    //--cros
    if(cmdOptionExists(argv, argv+argc, "--cros")){
        allow_cros = true;
        printf("CrownROS: %d \n", allow_cros);
    }
    
    //--centered-distribution
    if(cmdOptionExists(argv, argv+argc, "--centered-distribution")){
        centered_distribution = true;
        printf("Centered distribution (ellipse): %d \n", centered_distribution);
    }
    
    //--deactivate-offset
    if(cmdOptionExists(argv, argv+argc, "--deactivate-offset")){
        deactivate_offset = true;
        printf("Adjusted elliptical offsed deactivated: %d \n", deactivate_offset);
    }


    // Floats and ints
    // defaults
    int dsim_years = 1;
    int dnsims = 1;
    int dweather_period_len= 60;
    int dweather_files = 1;
    int dmax_fire_periods= 10000000;
    int dseed = 123;
    int diradius = 0;
    int delliptical_option=1;
    int dlb_formula=0;
    float dROS_Threshold= 0.0;
    float dHFI_Threshold= 0.0;
    float dCROS_Threshold= 0.5;
    float dCROSAct_Threshold= 1.0;
    float dROSCV= 0.;
    float dHFactor = 1.0;
    float dFFactor = 1.0;
    float dBFactor = 1.0;
    float dEFactor = 1.0;
    float dCBDFactor = 0.0;
    float dCBHFactor = 1.0;
    float dCCFFactor = 0.0;
    float dROS10Factor = 3.34;
    float dinput_PeriodLen= 1.;
    float dKFactor1= 1.;
    float dKFactor2= 2.;
    float dKFactor3= 0.;
    float dKFactor4= 0.;
    float dKFactor5= 0.;

    // aux 
    std::string::size_type sz;

    //--sim-years  (float)
    char * sim_years = getCmdOption(argv, argv + argc, "--sim-years");
    if (sim_years){
        printf("TotalYears: %s \n", sim_years);
        args_ptr->TotalYears = std::stoi (sim_years ,&sz); 
    }
    else args_ptr->TotalYears = dsim_years;

    //--nsims
    char * input_nsims = getCmdOption(argv, argv + argc, "--nsims");
    if (input_nsims){
        printf("TotalSims: %s \n", input_nsims);
        args_ptr->TotalSims = std::stoi (input_nsims ,&sz);  
    }
    else args_ptr->TotalSims = dnsims; 

    //--Weather-Period-Length
    char * weather_period_len = getCmdOption(argv, argv + argc, "--Weather-Period-Length");
    if (weather_period_len){
        printf("WeatherPeriodLength: %s \n", weather_period_len);
        args_ptr->MinutesPerWP = std::stoi (weather_period_len ,&sz); 
    }
    else args_ptr->MinutesPerWP = dweather_period_len;

    //--nweathers
    char * nweathers = getCmdOption(argv, argv + argc, "--nweathers");
    if (nweathers){
        printf("NWeatherFiles: %s \n", nweathers);
        args_ptr->NWeatherFiles = std::stoi (nweathers ,&sz); 
    }
    else args_ptr->NWeatherFiles = dweather_files;

    //--Fire-Period-Length
    char * input_PeriodLen = getCmdOption(argv, argv + argc, "--Fire-Period-Length");
    if (input_PeriodLen){
        printf("FirePeriodLength: %s \n", input_PeriodLen);
        if (std::stof (input_PeriodLen ,&sz) <= args_ptr->MinutesPerWP){
             args_ptr->FirePeriodLen = std::stof (input_PeriodLen ,&sz); 
        }
        else args_ptr->FirePeriodLen = args_ptr->MinutesPerWP;
    } 
    else args_ptr->FirePeriodLen = args_ptr->MinutesPerWP;

    //--IgnitionRad
    char * input_igrad = getCmdOption(argv, argv + argc, "--IgnitionRad");
    if (input_igrad){
        printf("IgnitionRadius: %s \n", input_igrad);
        args_ptr->IgnitionRadius = std::stoi (input_igrad ,&sz); 
    } 
    else args_ptr->IgnitionRadius =  diradius;
    
    //--EllipticalOption
    char * input_EllipticalOption = getCmdOption(argv, argv + argc, "--EllipticalOption");
    if (input_EllipticalOption){
        printf("EllipticalOption: %s \n", input_EllipticalOption);
        args_ptr->EllipticalOption = std::stoi (input_EllipticalOption ,&sz); 
    } 
    else args_ptr->EllipticalOption = delliptical_option;
    
    //--LBFormula
    char * input_LBFormula = getCmdOption(argv, argv + argc, "--LBFormula");
    if (input_LBFormula){
        printf("LBFormula: %s \n", input_LBFormula);
        args_ptr->LBFormula = std::stoi (input_LBFormula ,&sz); 
    } 
    else args_ptr->LBFormula = dlb_formula;

    
    //--ROS-Threshold
    char * ROS_Threshold = getCmdOption(argv, argv + argc, "--ROS-Threshold");
    if (ROS_Threshold){
        printf("ROSThreshold: %s \n", ROS_Threshold);
        args_ptr->ROSThreshold = std::stof (ROS_Threshold ,&sz); 
    }
    else args_ptr->ROSThreshold = dROS_Threshold;

    //--CROS-Threshold
    char * CROS_Threshold = getCmdOption(argv, argv + argc, "--CROS-Threshold");
    if (CROS_Threshold){
        printf("CROSThreshold: %s \n", CROS_Threshold);
        args_ptr->CROSThreshold = std::stof (CROS_Threshold ,&sz); 
    }
    else args_ptr->CROSThreshold = dCROS_Threshold;

    //--HFI-Threshold
    char * HFI_Threshold = getCmdOption(argv, argv + argc, "--HFI-Threshold");
    if (HFI_Threshold){
        printf("HFIThreshold: %s \n", HFI_Threshold);
        args_ptr->HFIThreshold = std::stof (HFI_Threshold ,&sz);
    }
    else args_ptr->HFIThreshold = dHFI_Threshold;

    //--CROSAct-Threshold
    char * CROSAct_Threshold = getCmdOption(argv, argv + argc, "--CROSAct-Threshold");
    if (CROSAct_Threshold){
        printf("CROSActThreshold: %s \n", CROSAct_Threshold);
        args_ptr->CROSActThreshold = std::stof (CROSAct_Threshold ,&sz); 
    }
    else args_ptr->CROSThreshold = dCROS_Threshold;


    //--HFactor
    char * H_Factor = getCmdOption(argv, argv + argc, "--HFactor");
    if (H_Factor){
        printf("HFactor: %s \n", H_Factor);
        args_ptr->HFactor = std::stof (H_Factor ,&sz); 
    }
    else args_ptr->HFactor = dHFactor;

    //---FFactor
    char * F_Factor = getCmdOption(argv, argv + argc, "--FFactor");
    if (F_Factor){
        printf("FFactor: %s \n", F_Factor);
        args_ptr->FFactor = std::stof (F_Factor ,&sz); 
    }
    else args_ptr->FFactor = dFFactor;

    ///--BFactor
    char * B_Factor = getCmdOption(argv, argv + argc, "--BFactor");
    if (B_Factor){
        printf("BFactor: %s \n", B_Factor);
        args_ptr->BFactor = std::stof (B_Factor ,&sz); 
    }
    else args_ptr->BFactor = dBFactor;

    ///--EFactor
    char * E_Factor = getCmdOption(argv, argv + argc, "--EFactor");
    if (E_Factor){
        printf("EFactor: %s \n", E_Factor);
        args_ptr->EFactor = std::stof (E_Factor ,&sz); 
    }
    else args_ptr->EFactor = dEFactor;

    ///--CBDFactor
    char * CBD_Factor = getCmdOption(argv, argv + argc, "--CBDFactor");
    if (CBD_Factor){
        printf("CBDFactor: %s \n", CBD_Factor);
        args_ptr->CBDFactor = std::stof (CBD_Factor ,&sz); 
    }
    else args_ptr->CBDFactor = dCBDFactor;

    ///--CCFFactor
    char * CCF_Factor = getCmdOption(argv, argv + argc, "--CCFFactor");
    if (CCF_Factor){
        printf("CCFFactor: %s \n", CCF_Factor);
        args_ptr->CCFFactor = std::stof (CCF_Factor ,&sz); 
    }
    else args_ptr->CCFFactor = dCCFFactor;

    ///--ROS10Factor
    char * ROS10_Factor = getCmdOption(argv, argv + argc, "--ROS10Factor");
    if (ROS10_Factor){
        printf("ROS10Factor: %s \n", ROS10_Factor);
        args_ptr->ROS10Factor = std::stof (ROS10_Factor, &sz); 
    }
    else args_ptr->ROS10Factor = dROS10Factor;

    //--ROS-CV
    char * ROS_CV = getCmdOption(argv, argv + argc, "--ROS-CV");
    if (ROS_CV){
        printf("ROS-CV: %s \n", ROS_CV);
        args_ptr->ROSCV = std::stof (ROS_CV ,&sz); 
    }
    else args_ptr->ROSCV = dROSCV;

    //--max-fire-periods
    char * max_fire_periods = getCmdOption(argv, argv + argc, "--max-fire-periods");
    if (max_fire_periods){
        printf("MaxFirePeriods: %s \n", max_fire_periods);
        args_ptr->MaxFirePeriods = std::stoi (max_fire_periods ,&sz); 
    }
    else args_ptr->MaxFirePeriods = dmax_fire_periods;

    //--seed  (int)
    char * seed = getCmdOption(argv, argv + argc, "--seed");
    if (seed){
        printf("seed: %s \n", seed);
        args_ptr->seed = std::stoi (seed ,&sz); 
    }
    else args_ptr->seed = dseed;
    
    // Elliptical factors to avoid long-shaped fires due to CA distortion
    //--KFactor1
    char * K_Factor_1 = getCmdOption(argv, argv + argc, "--KFactor1");
    if (K_Factor_1){
        printf("KFactor1: %s \n", K_Factor_1);
        args_ptr->KFactor1 = std::stof (K_Factor_1 ,&sz); 
    }
    else args_ptr->KFactor1 = dKFactor1;
    
    //--KFactor2
    char * K_Factor_2 = getCmdOption(argv, argv + argc, "--KFactor2");
    if (K_Factor_2){
        printf("KFactor2: %s \n", K_Factor_2);
        args_ptr->KFactor2 = std::stof (K_Factor_2 ,&sz); 
    }
    else args_ptr->KFactor2 = dKFactor2;
    
    //--KFactor3
    char * K_Factor_3 = getCmdOption(argv, argv + argc, "--KFactor3");
    if (K_Factor_3){
        printf("KFactor3: %s \n", K_Factor_3);
        args_ptr->KFactor3 = std::stof (K_Factor_3 ,&sz); 
    }
    else args_ptr->KFactor3 = dKFactor3;
    
    //--KFactor4
    char * K_Factor_4 = getCmdOption(argv, argv + argc, "--KFactor4");
    if (K_Factor_4){
        printf("KFactor4: %s \n", K_Factor_4);
        args_ptr->KFactor4 = std::stof (K_Factor_4 ,&sz); 
    }
    else args_ptr->KFactor4 = dKFactor4;
    
    //--KFactor5
    char * K_Factor_5 = getCmdOption(argv, argv + argc, "--KFactor5");
    if (K_Factor_5){
        printf("KFactor5: %s \n", K_Factor_5);
        args_ptr->KFactor5 = std::stof (K_Factor_5 ,&sz); 
    }
    else args_ptr->KFactor5 = dKFactor5;
    

    // Populate structure
    // Strings 
    if (input_folder == &empty){
        args_ptr->InFolder = ""; 
    }
    else args_ptr->InFolder = input_folder; 

    if (output_folder == &empty && input_folder != &empty){
        args_ptr->OutFolder = args_ptr->InFolder + "simOuts";
    } else if(output_folder == &empty && input_folder == &empty){
        args_ptr->OutFolder = "simOuts";
    } else if(output_folder != &empty && input_folder == &empty){
        args_ptr->OutFolder = output_folder;
    } else if(output_folder != &empty && input_folder != &empty){
        args_ptr->OutFolder = output_folder;
    }

    if (input_weather == &empty){
        args_ptr->WeatherOpt = "rows";
    }
    else{
        args_ptr->WeatherOpt = input_weather;
    }

    if (input_hplan == &empty){
        args_ptr->HarvestPlan = ""; 
    }
    else args_ptr->HarvestPlan = input_hplan; 

    // Optimal K factors from CSV (MHK)
    if (kopt == &empty){
        args_ptr->KOption = "";
    }
    else{
        args_ptr->KOption = kopt;
    }    

    // booleans
    args_ptr->OutMessages = out_messages;
    args_ptr->Trajectories = out_trajectories; 
    args_ptr->NoOutput = no_output;
    args_ptr->verbose = verbose_input; 
    args_ptr->Ignitions = input_ignitions;
    args_ptr->OutputGrids = out_grids; 
    args_ptr->FinalGrid = out_finalgrid;
    args_ptr->PromTuned = prom_tuned;
    args_ptr->Stats = out_stats;   
    args_ptr->BBOTuning = bbo_tuning;
    args_ptr->AllowCROS = allow_cros;
    args_ptr->CenteredDistribution = centered_distribution;
    args_ptr->NoAngleOffset = deactivate_offset;
    

}

// Print args
void printArgs(arguments args){
    std::cout << "InFolder: " << args.InFolder << std::endl; 
    std::cout << "OutFolder: " << args.OutFolder << std::endl; 
    std::cout << "WeatherOpt: " << args.WeatherOpt << std::endl;
    std::cout << "NWeatherFiles: " << args.NWeatherFiles << std::endl;	
    std::cout << "MinutesPerWP: " << args.MinutesPerWP << std::endl; 
    std::cout << "MaxFirePeriods: " << args.MaxFirePeriods << std::endl; 
    std::cout << "Messages: " << args.OutMessages << std::endl; 
    std::cout << "HarvestPlan: " << args.HarvestPlan << std::endl; 
    std::cout << "TotalYears: " << args.TotalYears << std::endl; 
    std::cout << "TotalSims: " << args.TotalSims << std::endl; 
    std::cout << "FirePeriodLen: " << args.FirePeriodLen << std::endl; 
    std::cout << "Ignitions: " << args.Ignitions << std::endl; 
    std::cout << "IgnitionRad: " << args.IgnitionRadius << std::endl; 
    std::cout << "OutputGrid: " << args.OutputGrids << std::endl; 
    std::cout << "FinalGrid: " << args.FinalGrid << std::endl; 
    std::cout << "PromTuned: " << args.PromTuned << std::endl; 
    std::cout << "BBOTuning: " << args.BBOTuning << std::endl; 
    std::cout << "EllipticalOption: " << args.EllipticalOption << std::endl; 
    std::cout << "CenteredDistribution: " << args.CenteredDistribution << std::endl; 
    std::cout << "LBFormula: " << args.LBFormula << std::endl; 
    std::cout << "NoAngleOffset: " << args.NoAngleOffset << std::endl; 
    std::cout << "Statistics: " << args.Stats << std::endl; 
    std::cout << "noOutput: " << args.NoOutput << std::endl; 
    std::cout << "verbose: " << args.verbose << std::endl; 
    std::cout << "seed: " << args.seed << std::endl; 



}