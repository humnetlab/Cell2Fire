// Include classes
#include "CellsFBP.h"
#include "SpottingFBP.h"
#include "fuelmodelBP.h"
#include "ReadCSV.h"
#include "ReadArgs.h"
#include "Ellipse.h"

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

// Torch
#include <torch/torch.h>
#include <torch/script.h> 

using namespace std;


/*************************
       Constructor   
*************************/
CellsFBP::CellsFBP(int _id, 
                   double _area, 
                   std::vector<int> _coord,  
                   int _fType, 
                   std::string _fType2,
                   double _perimeter, 
                   int _status,
                   std::unordered_map<std::string, int> & _adjacents, 
                   int _realId)
{
    // Global "dictionaries" (vectors) for status and types
    // Status: 0: "Available", 1: "Burning", 2: "Burnt", 3: "Harvested", 4:"Non Fuel"
    this->StatusD[0] = "Available";
    this->StatusD[1] = "Burning";
    this->StatusD[2] = "Burnt";
    this->StatusD[3] = "Harvested";
    this->StatusD[4] = "Non Fuel";

    // FTypeD: 0: "NonBurnable", 1: "Normal", 2: "Burnable
    this->FTypeD[0] = "NonBurnable";
    this->FTypeD[1] = "Normal";
    this->FTypeD[2] = "Burnable";

    // Initialize fields of the Cell object 
    this->id = _id;
    this->area = _area;
    this->coord = _coord;
    this->fType = _fType;
    this->fType2 = _fType2;
    this->perimeter = _perimeter;
    this->status = _status;
    this->adjacents = _adjacents;
    this->realId = _realId;
    this->_ctr2ctrdist = std::sqrt(this->area);

    if (std::abs(4 * this->_ctr2ctrdist - this->perimeter) > 0.01 * this->perimeter) {
        std::cerr << "Cell ID=" << this->id << "Area=" <<  this->area <<  "Perimeter=" <<  this->perimeter << std::endl;
        // maybe raise runtime error
    }
        
    // Inner fields
    this->gMsgList = std::unordered_map<int, std::vector<int>>();
    this->hPeriod = 0;

    this->fireStarts = 0;
    this->harvestStarts = 0;
    this->fireStartsSeason = 0;
    this->tYears = 4;

    this->gMsgListSeason = std::unordered_map<int, std::vector<int>>();
    this->fireProgress = std::unordered_map<int, double>();
    this->angleDict = std::unordered_map<int, double>();
    this->ROSAngleDir = std::unordered_map<int, double>();
    this->distToCenter = std::unordered_map<int, double>();
    this->angleToNb = std::unordered_map<int, int>();
}


/*******************************************************************
    Populates angles, distances, and initialize ROS per axis
    Modified by dlw to use cell area to compute distances in meters.
    ASSUME square fire cells.
    
    Returns        void
    
    Inputs:
    CoorCells      array of 2D int arrays
    AvailSet       int set
********************************************************************/
void CellsFBP::initializeFireFields(std::vector<std::vector<int>> & coordCells,    
                                    std::unordered_set<int> & availSet)
{  
    for (auto & nb : this->adjacents) {
        // CP Default value is replaced: None = -1
        if (nb.second != -1) {
            int a = -1 * coordCells[nb.second - 1][0] + coordCells[this->id][0];
            int b = -1 * coordCells[nb.second - 1][1] + coordCells[this->id][1];
            
            int angle = -1;
            if (a == 0) {
                if (b >= 0) 
                    angle = 270; 
                else 
                    angle = 90;
            }
            else if (b == 0) {
                if (a >= 0)
                    angle = 180;
                else
                    angle = 0;
            }
            else {
                // TODO: check this logi
                double radToDeg = 180 / M_PI;
                // TODO: i think all the negatives and abs cancel out
                double temp = std::atan(b * 1.0 / a) * radToDeg;
                if (a > 0)
                    temp += 180;
                if (a < 0 && b > 0) 
                    temp += 360;
                angle = temp;
            }

            this->angleDict[nb.second] = angle;
            if (availSet.find(nb.second) != availSet.end()) {
                this->ROSAngleDir[angle] = -1;
            }
            this->angleToNb[angle] = nb.second;
            this->fireProgress[nb.second] = 0.0;
            this->distToCenter[nb.second] = std::sqrt(a * a + b * b) * this->_ctr2ctrdist;
        }
    }
}


/******************************************************************************************
    Simple ROS allocation function for quick tests interpolating pairs of ROS and angles
    Returns      double

    Inputs:
    offset       double
    base         double
    ros1         double
    ros2         double
*******************************************************************************************/
double CellsFBP::allocate(double offset,
                          double base,
                          double ros1,
                          double ros2) {   
    
    // Simple interpolation
    double d = (offset - base) / 90;
    return (1 - d) * ros1 + d * ros2;
}
        
            
/*************************************************************************************
    New functions for calculating the ROS based on the fire angles
    Distribute the rate of spread (ROS,ros) to the axes given in the AngleList.
    All angles are w.t.r. E-W with East positive and in non-negative degrees.
    Inputs:
            thetafire: direction of "forward"
            forward : forward ROS
            flank: ROS normal to forward (on both sides)
            back: ROS in the opposide direction of forward
            AngleList: List of angles for the axes connecting centers
                       of interest (might have less than 8 angles)
    Effect:
            Populate the ROSAngleDir, whose indexes are the angles,
            with ROS values.
        
    Returns       void
*************************************************************************************/
void CellsFBP::ros_distr_old(double thetafire, 
                             double forward,
                             double flank,
                             double back) { 
    
    // Iterate through angles
    for (auto & angle : this->ROSAngleDir) {
        double offset = std::abs(angle.first - thetafire);
        
        double base = ((int)(offset)) / 90 * 90;
        double result;
        
        // Distribute ROS
        if (offset >= 0 && offset <= 90) {
            result = this->allocate(offset, 0, forward, flank);
        } else if (offset > 90 && offset < 180) {
            result = this->allocate(offset, 90, flank, back);
        } else if (offset >= 180 && offset <= 270) {
            result = this->allocate(offset, 180, back, flank);
        } else if (offset > 270 && offset < 360) {
            result = this->allocate(offset, 270, flank, forward);
        }
        this->ROSAngleDir[angle.first] = result;
    }
}



/********************************************************************
     Auxiliary function to round up numbers up to N decimals
********************************************************************/
double CellsFBP::round_up(double value, 
                          int decimal_places){
    
    // Constant
    const double multiplier = std::pow(10.0, decimal_places);
    
    // Return value
    return std::ceil(value * multiplier) / multiplier;
}



/************************************************************************************************************ 
    Estimates LB for different spread models given ws in [km/hr] and int indicating the scheme to use
    # 1. FBP System (Others)  
    # 2. FBP System (Grass)   
    # 3. Anderson (1983) - dense forest stand  
    # 4. Anderson (1983) - open forest stand   
    # 5. Anderson (1983) - grass/slash   
    # 6. Anderson (1983) - heavy slash   
    # 7. Alexander et al. (1985)   
    # 8. KITRAL System  
************************************************************************************************************/
double CellsFBP::general_lb(double ws_kmhr,
                            int lb_estimation_scheme){
    
    // Init
    double lb;
    double l_1, l_2;

    // Constants
    if (lb_estimation_scheme == 1){
        l_1 = 3.053;
        l_2 = 0.02667;
    } 
    else if (lb_estimation_scheme == 2){
        l_1 = 2.454;
        l_2 = 0.07154;
    }
    else if (lb_estimation_scheme == 3){
        l_1 = 1.411;
        l_2 = 0.01745;
    }
    else if (lb_estimation_scheme == 4){
        l_1 = 2.587;
        l_2 = 0.01142;
    }
    else if (lb_estimation_scheme == 5){
        l_1 = 5.578;
        l_2 = 0.0060230;
    }
    else if (lb_estimation_scheme == 6){
        l_1 = 3.749;
        l_2 = 0.0009885;
    }
    else if (lb_estimation_scheme == 7){
        l_1 = 3.063;
        l_2 = -0.01165;
    }
    else {
        l_1 = 2.233;
        l_2 = -0.01031;
    }
    
    // Formula
    lb = 1.0 + std::pow( l_1 * (1 - std::exp(-l_2 * (ws_kmhr*1.609))) , 2); // mph to kmhr conversion
    
    // Return estimated LB
    return lb;
}



/***************************************************************************************
    Parameterize ellipse from foci and calculates the expected ROS given an angle
    Inputs:
            thetafire: direction of "forward" (360 degrees)
            a: semi-major axis of the ellipse
            b: semi-minor axis of the ellipse
    Effect:
            Calculates ROS according to the ellipse 
        
    Returns:
            r: rho vector [ROS magnitude given theta]
    
****************************************************************************************/
double CellsFBP::rhoTheta(double theta, 
                          double a,
                          double b){
    
    // PI constant to get radians
    const double pi = 3.141592653589793;

    // Init
    double c2, e, r1, r2, r;

    // Calculations for RHO
    c2 = std::pow(a, 2) - std::pow(b, 2);
    e = std::sqrt(c2) / a;

    r1 = a * (1 - std::pow(e, 2));
    r2 = 1 - e * std::cos(theta * pi / 180.0);
    
    // RHO
    r = r1 / r2;

    // Return estimated RHO from foci as center
    return r;
}


/****************************************************************************************
    Parameterize ellipse from center and calculates the expected ROS given an angle
    Inputs:
            thetafire: direction of "forward" (360 degrees)
            a: semi-major axis of the ellipse
            b: semi-minor axis of the ellipse
    Effect:
            Calculates ROS according to the ellipse 
        
    Returns:
            r: rho vector [ROS magnitude given theta]
    
 ****************************************************************************************/
double CellsFBP::rhoTheta_center(double theta, 
                                 double a,
                                 double b){
    
    // PI constant to get radians
    const double pi = 3.141592653589793;
    
    // Init
    double theta_rads, r;
    
    // Radians
    theta_rads = theta * pi / 180.0;

    // Rho from center of ellipse
    r = (a * b) / std::sqrt(  (b * std::pow( std::sin(theta_rads), 2) ) + (a * std::pow( std::cos(theta_rads), 2) )  );

    // Return estimated ROS 
    return r;
}

/*********************************************************************************************************************
    Parameterize ellipse from foci including adjustment facotrs kopts to control elliptical shape
    Returns ROS given angle and previously calculated elliptical parameters
    Inputs:
            angle: angle to which we want to estimate ROS
            thetafire: direction of "forward" (360 degrees)
            LW: Length to bread ratio (using LW to follow paper)
            Ebar | Ebar_c: Eccentricities without and with kopts adjustments
            hros_mmin_c: head ros in meters / min
            k_opts: vector of optimal adjustment parameters 
            deactivate_angle_offset: boolean indicating if we deactivate the angle adjustment for back and front
            verbose: information level boolean
    Effect:
            Calculates ROS according to the ellipse 
        
    Returns:
            r: rho vector [ROS magnitude given theta]
    
 *********************************************************************************************************************/
double CellsFBP::rhoTheta_adjusted(double angle,
                                   double LW,
                                   double Ebar,
                                   double Ebar_c,
                                   double hros_mmin_c,
                                   std::vector<double> k_opts,
                                   bool deactivate_angle_offset,
                                   bool verbose){
    
    // PI constant to get radians
    const double pi = 3.141592653589793;
    
    // Original angles container and reference
    double theta_original = angle;
    double theta = angle;

    // Negative angles adjustment
    if (theta_original > 180){
        theta -= 360;
    }

    // Radians (explicit calculation)
    double theta_rad = theta * pi / 180.;
    double theta_c_rad;

    // Adjustment 
    if (std::abs(theta_rad) < pi / 2.){ 
        theta_c_rad = std::atan( std::tan(theta_rad) / std::pow(LW, k_opts[2]) );
    }

    else if (std::abs(theta_rad) > pi / 2.){ 
        theta_c_rad = std::atan( std::tan(theta_rad) / std::pow(LW, k_opts[3]) );
    }

    else{
        theta_c_rad = theta_rad;
    }

    // Adjustment to obtain the angle we want
    if (deactivate_angle_offset==false){
        if (theta < -90){
            theta_c_rad = theta_c_rad - pi;
        }
        if (theta > 90){
            theta_c_rad = theta_c_rad + pi;
        }
    }

    // Theta c in degrees [for reference]
    double theta_c = theta_c_rad * 180. / pi;

    // ROS calculations
    double Rtheta_c1 = hros_mmin_c * (1 - Ebar_c) / (1 - Ebar_c * std::cos(theta_c_rad));
    double Rtheta_c2 = hros_mmin_c * ( (1 - Ebar_c) / (1 + Ebar_c)  - (1 - Ebar) / (1 + Ebar) ) ;
    double Rtheta_c3 = std::pow( std::abs(theta_c_rad) / pi, k_opts[4] );
    double Rtheta_c = Rtheta_c1 - Rtheta_c2 * Rtheta_c3;

    // Return corrected angle and ROS(theta_corrected)
    return Rtheta_c;
}

/********************************************************************************************
    Elliptical fittings
    Alexander elliptical definition using HROS m/min and WS in km/hr
    b: semi-major | a: semi-minor | c: ignition to center (centre to foci)
********************************************************************************************/
std::vector<double> CellsFBP::alexander_ellipse(double hros_mmin, 
                                                double ws_kmhr,
                                                int lb_estimation_scheme,
                                                bool verbose){
    
    // Container
    std::vector<double> params;

    // Wind km/hr to m/s
    // double factor = 1000. / 3600.;
    double factor = 0.44704; // mph to m/s factor
    double ws_ms = ws_kmhr * factor;
    
    // Info
    if (verbose){
        std::cout << "ws Km/h: " << ws_kmhr << std::endl;
        std::cout << "ws m/s: " << ws_ms << std::endl;
        std::cout << "hros m/min: " << hros_mmin << std::endl;
    }

    // Fit ellipse
    double LB, HB, a, b, c;
    if (lb_estimation_scheme == 0){ // FARSITE
        LB = 0.936 * std::exp(0.2566 * ws_ms) + 0.461 * std::exp(-0.1548 * ws_ms) - 0.397;
    }
    else if (lb_estimation_scheme == -1){ // Behave
        LB = (0.936 * std::exp(0.2566 * ws_ms) + 0.461 * std::exp(-0.1548 * ws_ms) - 0.397) * std::exp(0.46);
    }
    else{
        LB = general_lb(ws_kmhr, lb_estimation_scheme);
    }
    
    HB = (LB + std::pow(std::pow(LB, 2) - 1, 0.5)) / (LB - std::pow(std::pow(LB, 2) - 1, 0.5));
    a = 0.5 * (hros_mmin + hros_mmin / HB) / LB;
    b = (hros_mmin + hros_mmin / HB) / 2;
    c = b - hros_mmin / HB;

    // Info
    if (verbose){
        std::cout << "Semi-major (b): " << b << std::endl;
        std::cout << "Semi-minor (a): " << a << std::endl;
        std::cout << "c: " << c << std::endl;
    }
        
    // Record
    params.push_back(a);
    params.push_back(b);
    params.push_back(c);
    
    // Return ellipse values
    //return a, b, c;
    return params;
}

/********************************************************************************************
           FBP ellipse using HROS, FROS and BROS in the same unit [e.g., m/min]
********************************************************************************************/
std::vector<double> CellsFBP::fbp_ellipse(double forward,
                                          double flank,
                                          double back,
                                          bool verbose) {  
    
    // Container
    std::vector<double> params;
 
    // Initialize parameters
    double a = (forward + back) / 2;
    double b;
    
    // Create vectors for fitting
    std::vector<double> _x = {0.0, back, back, (forward + back) / 2., (forward + back) / 2., (forward + back)};
    std::vector<double> _y = {0.0, std::pow(flank, 2) / a, - (std::pow(flank, 2) / a), flank, -flank, 0.0};

    // Fit the Ellipse calling ellipse class
    Ellipse SqlEllipse(_x, _y);     
    
    // Get parameters
    params = SqlEllipse.get_parameters();
    a = params[0];
    b = params[1];
    
    // INFO
    if (verbose){
        std::cout << "a: " << a << std::endl; 
        std::cout << "b: " << b << std::endl; 
    }
 
    // Return ellipse values (c not needed)
    return params;
}

/********************************************************************************************
    Basic elliptical distribution using provided HROS m/min and LB
    b: semi-major | a: semi-minor | c: ignition to center (centre to foci)
********************************************************************************************/
std::vector<double> CellsFBP::lb_ellipse(double ros_mmin, 
                                         double LB,
                                         bool verbose){
    
    // Container
    std::vector<double> params;
 
    // Info
    if (verbose){
        std::cout << "lb: " <<  LB << std::endl;
        std::cout << "ros [m/min]: " << ros_mmin << std::endl;
    }
    
    // Init
    double HB, a, b, c;
    
    // Calculate
    HB = (LB + std::pow(std::pow(LB, 2) - 1, 0.5)) / (LB - std::pow(std::pow(LB, 2) - 1, 0.5));
    a = 0.5 * (ros_mmin + ros_mmin / HB) / LB;
    b = (ros_mmin + ros_mmin / HB) / 2;
    c = b - ros_mmin / HB;
        
    // Info
    if (verbose){
        std::cout << "Semi-major (b): " << b << std::endl;
        std::cout << "Semi-minor (a): " << a << std::endl;
        std::cout << "c: " << c << std::endl;
    }
     
    
    // Record
    params.push_back(a);
    params.push_back(b);
    params.push_back(c);
    
    // Return ellipse values
    //return a, b, c;
    return params;
}
    

/********************************************************************************************
    Estimate elliptical parameters as a function of the fire main direction/angle 
    HROS in [m/min], ws in [km/hr], and adjusting factors for elliptical propagation
********************************************************************************************/
std::vector<double> CellsFBP::ellipse_adjusted(double hros_mmin,
                                               double ws_kmhr,
                                               std::vector<double> k_opts, // unordered map, key=float, value=std::vector
                                               std::unordered_map<double, std::vector<double>> *kdict_ptr,
                                               int lb_estimation_scheme,
                                            //    arguments * args,
                                               std::string KOption,
                                               bool verbose){
    
    // Container
    std::vector<double> params;
    
    // Initialize
    double ws_ms, LW, Ebar_c, Ebar, ros_mmin_c, a, b, c;
    // std::string KOption;

    // Wind speed transformation from km/hr to m/s
    // ws_ms = (ws_kmhr * 1000.) / 3600.;
    ws_ms = ws_kmhr * 0.44704; // mph to m/s
    
    // std::cout << "DEBUGGING ellipse_adjused : " << lb_estimation_scheme << std::endl;

    // LB = LW using alpha and beta parameters
    if (lb_estimation_scheme == 0){
        LW = 0.936 * std::exp(0.2566 * ws_ms) + 0.461 * std::exp(-0.1548 * ws_ms) - 0.397;
    }
    else{
        LW = general_lb(ws_kmhr, lb_estimation_scheme);
    }
    
    // Eccentricities
    Ebar = std::sqrt( 1. - 1. / std::pow(LW, 2) );

    // Pass the library (table) and resample from k opt library
    // int numRows = sizeof(kdict_ptr) / sizeof(kdict_ptr[0]); // get the number of rows
    // int numCols = sizeof(kdict_ptr[0]) / sizeof(double); // get the number of columns

    // std::cout << "DEBUGGING Ellipse_adjusted" << std::endl;
    // initialize for args flag
    std::string Empty = "";
    const char * EM = Empty.c_str(); 

    // std::unordered_map<double, std::vector<double>> *kdict_ptr = new std::unordered_map<double, std::vector<double>>;
    // std::cout << "Debugging in CellsFBP ellipse adjusted: Check to see contents of kdict_ptr" << std::endl;
    // for (auto iter = kdict_ptr->begin(); iter != kdict_ptr->end(); ++iter) {
    //     std::cout << iter->first << ": [";
    //     for (double d : iter->second) {
    //         std::cout << d << ", ";
    //     }
    //     std::cout << "]\n";
    // }

    // std::cout << "##### INPUT K FACTORS BEFORE LOOP #####" << std::endl;    
    // std::cout << "K1 = " << k_opts[0] << std::endl;
    // std::cout << "K2 = " << k_opts[1] << std::endl;
    // std::cout << "K3 = " << k_opts[2] << std::endl;
    // std::cout << "K4 = " << k_opts[3] << std::endl;
    // std::cout << "K5 = " << k_opts[4] << std::endl;
    
    // std::cout << "KOPTION : " << KOption << std::endl;

    if(strcmp(KOption.c_str(), EM) != 0){
    // if (KOption=='csv'){
        // std::cout << "DEBUGGING Ellipse_adjusted with K Opt CSV" << std::endl;
        std::vector<double> k_opts;
        // std::unordered_map<double, std::vector<double>> koptDict;
        // kdict = &kdict_ptr

        double closestIndex = 0; // initialize the closest index variable (Index = Eccentricity)
        // std::cout << "DEBUGGING Ellipse_adjusted Entered loop 0_2" << std::endl;
        
        if (!kdict_ptr) {
            std::cerr << "Error: empty pointer." << std::endl;
        }

        // for (auto iter = kdict_ptr->begin(); iter != kdict_ptr->end(); ++iter) {
        //     std::cout << iter->first << ": [";
        //     for (double d : iter->second) {
        //         std::cout << d << ", ";
        //     }
        //     std::cout << "]\n";
        // }

        double closestValue = kdict_ptr->begin()->first; // initialize to the first key in the dictionary
        double closestDistance = std::abs(closestValue - Ebar); // initialize to the absolute difference between the first value and the input

        for (auto it = kdict_ptr->begin(); it != kdict_ptr->end(); ++it) {
            double currentDistance = std::abs(it->first - Ebar);
            if (currentDistance < closestDistance) {
                closestIndex = it->first;
                closestDistance = currentDistance;
            }
        }

        // Dictionary order: F1-score, K1, K2, K3, K4, K5
        // std::cout << "Size of K opts vector : " << k_opts.size() << std::endl;
        k_opts.push_back((*kdict_ptr)[closestIndex][1]);
        k_opts.push_back((*kdict_ptr)[closestIndex][2]);
        k_opts.push_back((*kdict_ptr)[closestIndex][3]);
        k_opts.push_back((*kdict_ptr)[closestIndex][4]);
        k_opts.push_back((*kdict_ptr)[closestIndex][5]);
        // std::cout << "DEBUGGING Ellipse_adjusted Kopts vector done" << std::endl;
        
        if (verbose){
        std::cout << "##### OPTIMIZED K FACTORS #####" << std::endl;
        std::cout << "K1 = " << k_opts[0] << std::endl;
        std::cout << "K2 = " << k_opts[1] << std::endl;
        std::cout << "K3 = " << k_opts[2] << std::endl;
        std::cout << "K4 = " << k_opts[3] << std::endl;
        std::cout << "K5 = " << k_opts[4] << std::endl;
        }
        // std::cout << "CellFBP error : DEBUGGING check" << std::endl;
        
        Ebar_c = std::sqrt( 1. - 1. / std::pow(LW, k_opts[1]) );

        // std::cout << "CellFBP error : DEBUGGING check2" << std::endl;
        // Adjusted HROS
        ros_mmin_c = hros_mmin * k_opts[0];
        // std::cout << "CellFBP error : DEBUGGING check3" << std::endl;  
        params.push_back(a);
        params.push_back(b);
        params.push_back(c);
        params.push_back(LW);
        params.push_back(Ebar);
        params.push_back(Ebar_c);
        params.push_back(ros_mmin_c);

        params.push_back(k_opts[0]);
        params.push_back(k_opts[1]);
        params.push_back(k_opts[2]);
        params.push_back(k_opts[3]);
        params.push_back(k_opts[4]);
    }
    else {

    // std::cout << "DEBUGGING Ellipse_adjusted WITHOUT K Opt CSV!!!" << std::endl;
    Ebar_c = std::sqrt( 1. - 1. / std::pow(LW, k_opts[1]) );
    ros_mmin_c = hros_mmin * k_opts[0];

    // Record
    // TODO : Pass k_opts to params
    params.push_back(a);
    params.push_back(b);
    params.push_back(c);
    params.push_back(LW);
    params.push_back(Ebar);
    params.push_back(Ebar_c);
    params.push_back(ros_mmin_c);

    // std::cout << "DEBUGGING OUT OF LOOP" << std::endl;
    // std::cout << "K1 = " << k_opts[0] << std::endl;
    // std::cout << "K2 = " << k_opts[1] << std::endl;
    // std::cout << "K3 = " << k_opts[2] << std::endl;
    // std::cout << "K4 = " << k_opts[3] << std::endl;
    // std::cout << "K5 = " << k_opts[4] << std::endl;    

    }            
    // Info
    if (verbose){
    std::cout << "ws [Km/h]: " << ws_kmhr << std::endl;
    std::cout << "ws [m/s]: " << ws_ms << std::endl;
    std::cout << "ros [m/min]: " << hros_mmin << std::endl;
    std::cout << "ros_c [m/min]: " << round_up(ros_mmin_c, 5) << std::endl; 
    std::cout << "LB: " << LW << std::endl;
    std::cout << "E: " << Ebar << std::endl;
    std::cout << "E_c: " << Ebar_c << std::endl;
    }    

    // Return relevant estimations
    return params;
}


/********************************************************************************************
    Distributes ROS values (i.e., rho vector values = rate of spread) at different angles 
    given wind speed, direction, 
    Multiple inputs depending on the elliptical scheme selected
********************************************************************************************/
void CellsFBP::ros_distr(double thetafire,
                         double forward,
                         double flank,
                         double back,
                         double lb,
                         double ws,
                         int elliptical_option,
                         bool centered_distribution,
                         bool deactivate_angle_offset,
                         double EFactor,
                        //  std::unordered_map<double, std::vector<double>> koptDict,
                         std::unordered_map<double, std::vector<double>> *kdict_ptr,
                         std::vector<double> k_opts,
                         std::vector<double> ros_boosters,
                         int lb_estimation_scheme,
                         arguments * args,
                         bool verbose) {   
    
    // Containers
    double a, b, c;
    double Ebar, Ebar_c, ros_mmin_c;
    std::vector<double> params;

    // Args
    // this->args = &this->args;
    
    // Estimate a, b from ellipse
    // Alexander by default
    if (elliptical_option == 1){
        if (verbose){
            std::cout << "\nAlexander Distribution" << std::endl;
        }
        params = alexander_ellipse(forward, ws, lb_estimation_scheme, verbose);
        a = params[0];
        b = params[1];
        
    // FBP
    } else if (elliptical_option == 2){ 
        if (verbose){
            std::cout << "\nFBP Distribution" << std::endl;
        }
        params = fbp_ellipse(forward, flank, back, verbose);
        b = params[0];
        a = params[1];
        
    // LB provided (coming from the BP calculate function)
    } else if (elliptical_option == 3){
        if (verbose){
            std::cout << "\nLB Distribution" << std::endl;
        }
        params = lb_ellipse(forward, lb, verbose);
        a = params[0];
        b = params[1];
        
    // Adjusted Alexander    
    } else if (elliptical_option == 4){
        if (verbose){
            std::cout << "\nAdjusted Alexander Distribution" << std::endl;
        }

        // if (!kdict_ptr) {
        //     std::cerr << "ros_distr Error: The input pointer is null." << std::endl;
        //     return;
        // }

        // initialize for args flag
        // std::cout << "DEBUGGING in ros_distr" << std::endl;
        // std::cout << "CellsFBP KOPTION : " << args->KOption << std::endl;

        // std::string Empty = "";
        // const char * EM = Empty.c_str(); 

        // Doing the call for k_opts dictionary k_opt[Ebar] // MHK 021723
        // std::cout << "DEBUGGING in ros_distr before loop" << std::endl;
        
        if (args->KOption.empty()) {
        // if(strcmp(args->KOption.c_str(), EM) != 0){
            // std::vector<double> k_opts;
            // std::unordered_map<double, std::vector<double>> koptDict = new std::unordered_map<double, std::vector<double>>;

            // std::cout << "DEBUGGING in ros_distr ing" << std::endl;
            // std::unordered_map<double, std::vector<double>> koptDict;
            // std::unordered_map<double, std::vector<double>> *kdict_ptr = new std::unordered_map<double, std::vector<double>>;

            if (!kdict_ptr) {
                std::cerr << "before params Error: kdict_ptr is empty." << std::endl;
                return;
            }
            
            // std::cout << "Debugging in CellsFBP : Check to see contents of kdict_ptr" << std::endl;
            // for (auto iter = kdict_ptr->begin(); iter != kdict_ptr->end(); ++iter) {
            //     std::cout << iter->first << ": [";
            //     for (double d : iter->second) {
            //         std::cout << d << ", ";
            //     }
            //     std::cout << "]\n";
            // }

            params = ellipse_adjusted(forward, ws, k_opts, kdict_ptr, lb_estimation_scheme, args->KOption, verbose);            

        }
        
        params = ellipse_adjusted(forward, ws, k_opts, kdict_ptr, lb_estimation_scheme, args->KOption, verbose);            
        a = params[0];
        b = params[1];
        c = params[2]; 
        lb = params[3]; 
        Ebar = params[4]; 
        Ebar_c = params[5];
        ros_mmin_c = params[6];
    
        // TODO : Pass k opt vectors into params to read here
    

    // If no method, circle
    } else{
        // Circle if no option
        a = 1.0;
        b = 1.0;
        c = 0.0;
    }
    
    // Basic info from the 
    if (verbose){
        std::cout << "ROS DIST (a): " << a << std::endl;
        std::cout << "ROS DIST (b): " << b << std::endl;
    }
    
    // Ros allocation for each angle inside the dictionary
    for (auto & angle : this->ROSAngleDir) {
        double offset = angle.first - thetafire;

        // Adjust angle due to offset
        if (offset < 0) {
            offset += 360;
        }
        if (offset > 360) {
            offset -= 360;
        }

        // Centered or foci distributed
        if (centered_distribution){
            this->ROSAngleDir[angle.first] = rhoTheta_center(offset, a, b) * EFactor;
            
        // Foci
        } else{
            // Adjusted
            if (elliptical_option == 4){
                double _ros = rhoTheta_adjusted(offset,
                                                lb,
                                                Ebar,
                                                Ebar_c,
                                                ros_mmin_c,
                                                k_opts, // TODO : Should be reflecting the optimal corrected 
                                                deactivate_angle_offset,
                                                verbose);
                this->ROSAngleDir[angle.first] = _ros * EFactor;
            }
            
            else{
                this->ROSAngleDir[angle.first] = rhoTheta(offset, b, a) * EFactor;
            }
        }
    }

}



/********************************************************************************************
    Slope effect
    Inputs: 
        elev_i: elevation of burning cell
        elev_j: elevation of cell reached by fire
        cellsize: side of a cell 
/********************************************************************************************/

double CellsFBP::slope_effect(double elev_i, double elev_j, double cellsize)
  {
    double ps_ij = (elev_j - elev_i) / cellsize;
    double se;
    se = 1. + 0.023322 * ps_ij + 0.00013585 * std::pow(ps_ij, 2);
    
    return se;
  }



/********************************************************************************************
    Returns           vect[integers]   
    Important: we are sending a True sometimes, pick a special value -x for replacing it
    
    Inputs:
    period            int
    AvailSet          int set
    verbose           boolean
    df                Data frame
    coef              pointer
    spotting          boolean
    SpottingParams    Data frame
    CoordCells        array of 2D doubles arrays
    Cells_Obj         dictionary of cells objects
/********************************************************************************************/
std::vector<int> CellsFBP::manageFire(int period, std::unordered_set<int> & AvailSet,      
                                      inputs df_ptr[], 
                                      fuel_coefs * coef, 
                                    //   koptDF * kopt_ptr,
                                      std::unordered_map<double, std::vector<double>> *kdict_ptr,
                                      std::vector<std::vector<int>> & coordCells,
                                      std::unordered_map<int, CellsFBP> & Cells_Obj, 
                                      arguments * args,
                                      weatherDF * wdf_ptr, 
                                      std::vector<double> * FSCell,
                                      double randomROS) {
    // Special flag for repetition (False = -99 for the record)
    int repeat = -99;

    // msg lists contains integers (True = -100)
    std::vector<int> msg_list_aux;
    msg_list_aux.push_back(0);
    std::vector<int> msg_list;

    // Populate Inputs 
    df_ptr[this->realId-1].waz = wdf_ptr->waz;
    df_ptr[this->realId-1].ws = wdf_ptr->ws;
    df_ptr[this->realId-1].mc1 = wdf_ptr->mc1;
    df_ptr[this->realId-1].mc10 = wdf_ptr->mc10;
    df_ptr[this->realId-1].mc100 = wdf_ptr->mc100;
    df_ptr[this->realId-1].mcWoody = wdf_ptr->mcWoody;
    df_ptr[this->realId-1].mcHerb = wdf_ptr->mcHerb;
    df_ptr[this->realId-1].factor_cbd = args->CBDFactor;   
    df_ptr[this->realId-1].factor_ccf = args->CCFFactor;
    df_ptr[this->realId-1].factor_ros10 = args->ROS10Factor;
    df_ptr[this->realId-1].factor_actv = args->CROSActThreshold;
    df_ptr[this->realId-1].cros = args->AllowCROS;
    df_ptr[this->realId-1].verbose = args->verbose;

    // Compute main angle and ROSs: forward, flanks and back
    main_outs mainstruct;
    snd_outs sndstruct;
    fire_struc headstruct, backstruct, flankstruct;

    // Calculate parameters
    // calculate(&df_ptr[this->realId-1], coef, &mainstruct, &sndstruct, &headstruct, &flankstruct, &backstruct);
    calculate_BP(&df_ptr[this->realId-1], coef, &mainstruct, &sndstruct, &headstruct, &flankstruct, &backstruct);

    /*  ROSs DEBUG!   */
    if(args->verbose){
        std::cout << "*********** ROSs debug ************" << std::endl;
        std::cout <<  "-------Input Structure--------" << std::endl;
        std::cout <<  "fueltype: " << df_ptr[this->realId-1].fueltype << std::endl;
        std::cout <<  "nfueltype: " << df_ptr[this->realId-1].nftype << std::endl;
        std::cout <<  "ws: " << df_ptr[this->realId-1].ws << std::endl;
        std::cout <<  "waz: " << df_ptr[this->realId-1].waz << std::endl;
        std::cout <<  "ps: " << df_ptr[this->realId-1].ps << std::endl;
        std::cout <<  "saz: " << df_ptr[this->realId-1].saz << std::endl;
        std::cout <<  "cur: " << df_ptr[this->realId-1].cur << std::endl;
        std::cout <<  "elev: " << df_ptr[this->realId-1].elev << std::endl;
        std::cout <<  "cbd: " << df_ptr[this->realId-1].cbd << std::endl;
        std::cout <<  "cbh: " << df_ptr[this->realId-1].cbh << std::endl;
        std::cout <<  "ccf: " << df_ptr[this->realId-1].ccf << std::endl;
        std::cout <<  "cros: " << df_ptr[this->realId-1].cros << std::endl;
        std::cout << "HFactor: " << args->HFactor << std::endl;
        std::cout << "FFactor: " << args->FFactor << std::endl;
        std::cout << "Bactor: " << args->BFactor << std::endl;

        std::cout <<  "\n-------Mainout Structure--------" << std::endl;
        std::cout << "rss: " << mainstruct.rss << std::endl;
        std::cout << "angle: " << mainstruct.angle << std::endl;
        std::cout << "cros: " << mainstruct.cros << std::endl;
        std::cout <<  "\n-------Headout Structure--------" << std::endl;
        std::cout <<  "hros: " << headstruct.ros << std::endl;
        std::cout <<  "ros with Hfactor: " << headstruct.ros * args->HFactor << std::endl;
        std::cout <<  "rss: " << headstruct.rss << std::endl;
        std::cout <<  "\n------- Flank Structure--------" << std::endl;
        std::cout <<  "fros: " << flankstruct.ros << std::endl;
        std::cout <<  "ros with Ffactor: " << flankstruct.ros * args->FFactor << std::endl;
        std::cout <<  "rss: " << flankstruct.rss << std::endl;
        std::cout <<  "\n------- Back Structure--------" << std::endl;
        std::cout <<  "bros: " << backstruct.ros  << std::endl;
        std::cout <<  "ros with Bfactor " << backstruct.ros * args->BFactor << std::endl;
        std::cout <<  "rss: " << backstruct.rss << std::endl;
        std::cout <<  "\n------- Extra --------" << std::endl;
        std::cout <<  "lb: " << sndstruct.lb * args->BFactor << std::endl;
        std::cout << "a RSS: " << mainstruct.a << std::endl;
        std::cout << "b RSS: " << mainstruct.b << std::endl;
        std::cout << "c RSS: " << mainstruct.c << std::endl;
        std::cout << "*********** ROSs debug ************" << std::endl;
        }
    
    // Cartesian angle
    double cartesianAngle = wdf_ptr->waz;

    // Random variable for stochastic ROS
    double ROSRV = 0;
    if (args->ROSCV > 0) {
        ROSRV = randomROS;
    }

    // Display if verbose True (FBP ROSs, Main angle, and ROS std (if included))
    if (args->verbose) {
        std::cout << "Main Angle (raz): " << wdf_ptr->waz  << " Cartesian: " << cartesianAngle << std::endl;
        std::cout << "Front ROS Value: " << headstruct.rss * args->HFactor << std::endl; 
        std::cout << "Flanks ROS Value: " << flankstruct.rss * args->FFactor << std::endl;
        std::cout << "Rear ROS Value: " << backstruct.rss * args->BFactor << std::endl;
        std::cout << "Std Normal RV for Stochastic ROS CV: " << ROSRV << std::endl;
    }

    // If cell cannot send (thresholds), then it will be burned out in the main loop
    double HROS = (1. + args->ROSCV * ROSRV) * headstruct.rss * args->HFactor;
    
    // Extra debug step for sanity checks
    if (args->verbose){
        std::cout << "\nSending message conditions" << std::endl;
        std::cout << "HROS: " << HROS << " Threshold: "<<  args->ROSThreshold << std::endl;
    }

    // Check thresholds for sending messages
    if (HROS > args->ROSThreshold) {
        repeat = -100;
        
        if (args->verbose) {
            std::cout << "\nRepeat condition: " << repeat << std::endl;
            std::cout << "Cell can send messages" << std::endl;
        }
        
        // **ROS distribution method**
        // Get wind
        double wss = wdf_ptr->ws;
        

        // 1. Read and load K opt csv file + create dictionary (MHK)
        // Vector of k_opts
        std::vector<double> k_opts;  

        // if args exist:
        // Load K opt table and variables
        // double eccentricity = kopt_ptr.ecc
        // k_opts.push_back(kopt_ptr.K1);
        // k_opts.push_back(kopt_ptr.K2);
        // k_opts.push_back(kopt_ptr.K3);
        // k_opts.push_back(kopt_ptr.K4);
        // k_opts.push_back(kopt_ptr.K5);

        // else:

        // std::cout << "CellsFBP DEBUGGING Push back K opts" << std::endl;
        // MHK TO DO: Add flag if KFactors 1-5 exist and create default values if not existing
        if (args->KFactor1) {
            k_opts.push_back(args->KFactor1);
            k_opts.push_back(args->KFactor2);
            k_opts.push_back(args->KFactor3);
            k_opts.push_back(args->KFactor4);
            k_opts.push_back(args->KFactor5);
        }


        // ROS Boosters (placeholder in case they are needed in the future, for each neighbor)
        std::vector<double> ros_boosters;
        ros_boosters.push_back(1.);
        ros_boosters.push_back(1.);
        ros_boosters.push_back(1.);
        ros_boosters.push_back(1.);
        ros_boosters.push_back(1.);
        ros_boosters.push_back(1.);
        ros_boosters.push_back(1.);
        ros_boosters.push_back(1.);
        
        // std::cout << "Debugging in CellsFBP : Check to see contents of kdict_ptr" << std::endl;
        // for (auto iter = kdict_ptr->begin(); iter != kdict_ptr->end(); ++iter) {
        //     std::cout << iter->first << ": [";
        //     for (double d : iter->second) {
        //         std::cout << d << ", ";
        //     }
        //     std::cout << "]\n";
        // }
        
        // std::cout << "CellsFBP DEBUGGING Declare new var to derefernce pointer" << std::endl;
        // // Declare variable obtain dereferenced pointer (containing dictionary of K opt values)
        // // if (koptDict.empty()) {
        // if (koptDict.empty()) {
        //     std::cerr << "CellsFBP Error: The input pointer is null." << std::endl;
        //     // return;
        // }
        // if (!kdict_ptr) {
        //     std::cerr << "CellsFBP Error: The input pointer !kdict_ptr is null." << std::endl;
        //     // return;
        // }
        
        // std::unordered_map<double, std::vector<double>> *kdict_ptr = new std::unordered_map<double, std::vector<double>>;


        // if (koptDict.empty()) {
        //     std::cerr << "CellsFBP Error after: The input pointer is null." << std::endl;
        //     // return;
        // }
        // if (!kdict_ptr) {
        //     std::cerr << "CellsFBP Error after: The input pointer !kdict_ptr is null." << std::endl;
        //     // return;
        // }

        // std::cout << "CellsFBP DEBUGGING ros_distr" << std::endl;
        // Distribute ROS per theta angle
        ros_distr(cartesianAngle,
                  headstruct.ros * args->HFactor,
                  flankstruct.ros * args->FFactor, 
                  backstruct.ros * args->BFactor,
                  sndstruct.lb,
                  wss,
                  args->EllipticalOption,
                  args->CenteredDistribution,
                  args->NoAngleOffset,
                  args->EFactor,
                  kdict_ptr, // From MangeFire --> ros_distr --> ellipse_adjusted (MHK)
                  k_opts, // Table (Dataformat) format coming from main
                  ros_boosters,
                  args->LBFormula,
                  args,
                //   args->KOption,
                  args->verbose);
        
        // std::cout << "CellsFBP DEBUGGING ros_distr finished" << std::endl;
        
        // this is a iterator through the keyset of a dictionary
        for (auto&  _angle : this->ROSAngleDir) {
            double angle = _angle.first;
            int nb = angleToNb[angle];
            double ros = (1. + args->ROSCV * ROSRV) * _angle.second;
            
            // Clean nan ROS if numerical issues
            if(std::isnan(ros)){
                ros = 1e-4;
            }

            if (args->verbose) {
                std::cout << "     (angle, realized ros in m/min): (" << angle << ", " << ros << ")" << std::endl;
            }

            // Slope effect 
            double se = slope_effect(df_ptr[this->realId - 1].elev, df_ptr[nb-1].elev,  this->perimeter / 4.);
            if (args->verbose) { 
                std::cout << "Slope effect: " << se << std::endl;
            }

            // PeriodLen minutes and ros m/min
            this->fireProgress[nb] += ros * args->FirePeriodLen * se;   // Updates fire progress
            if (args->verbose) {
                std::cout << "Fire Progress : " << this->fireProgress[nb] << std::endl;
                std::cout << "Distance to Center : " << this->distToCenter[nb] << std::endl;
            }

            // If the message arrives to the adjacent cell's center, send a message
            if (this->fireProgress[nb] >= this->distToCenter[nb]) {
                msg_list.push_back(nb);
                FSCell->push_back(double(this->realId));
                FSCell->push_back(double(nb));
                FSCell->push_back(double(period));
                FSCell->push_back(ros);
                FSCell->push_back(angle);

            }    
            if (args->verbose) {
                std::cout << "DEBUGGING repeat : " << repeat << std::endl;
                std::cout << "DEBUGGING msg_list_aux : " << msg_list_aux[0] << std::endl;
            }

            // Info for debugging status of the cell and fire evolution
            if (this->fireProgress[nb] < this->distToCenter[nb] && repeat == -100 && -100 != msg_list_aux[0]){
                    if (args->verbose){
                        std::cout << "A Repeat = TRUE flag is sent in order to continue with the current fire....." << std::endl;
                        std::cout << "Main workaround of the new sim logic....." << std::endl;
                    }
                    msg_list_aux[0] = repeat;
            }

        }
    
        // Info fire progress
        if (args->verbose){
            for (auto & nb : this->fireProgress){
                std::cout << " " << nb.first << " : " << nb.second;
            }
            std::cout << std::endl;
        }
    }

    // If original is empty (no messages but fire is alive if aux_list is not empty)
    if  (msg_list.size() == 0){
        if (msg_list_aux[0] == -100){
            msg_list = msg_list_aux;
        }
        else{
            this->status = 2;   // we are done sending messages, call us burned
        }
    }
    
    // Info
    if (args->verbose){
        std::cout << " ----------------- End of new manageFire function -----------------" << std::endl;
    }
    
    // Messages list
    return msg_list;
}
    




/*********************************************************

    Manage fire for BBO tuning version (Not up to date)

**********************************************************/
std::vector<int> CellsFBP::manageFireBBO(int period, std::unordered_set<int> & AvailSet,      
                                         inputs * df_ptr,
                                         fuel_coefs * coef, 
                                        //  koptDF * kopt_ptr,
                                         std::unordered_map<double, std::vector<double>> *kdict_ptr,
                                         std::vector<std::vector<int>> & coordCells,
                                         std::unordered_map<int, CellsFBP> & Cells_Obj, 
                                         arguments * args,
                                         weatherDF * wdf_ptr, 
                                         std::vector<double> * FSCell,
                                         double randomROS,
                                         std::vector<float> & EllipseFactors) {
    // Special flag for repetition (False = -99 for the record)
    int repeat = -99;

    // msg lists contains integers (True = -100)
    std::vector<int> msg_list_aux;
    msg_list_aux.push_back(0);
    std::vector<int> msg_list;

    // Compute main angle and ROSs: forward, flanks and back
    main_outs mainstruct;
    snd_outs sndstruct;
    fire_struc headstruct, backstruct, flankstruct;

    // Populate inputs 
    df_ptr->waz = wdf_ptr->waz;
    df_ptr->ws = wdf_ptr->ws;
    df_ptr->factor_cbd = args->CBDFactor;   
    df_ptr->factor_ccf = args->CCFFactor;
    df_ptr->factor_ros10 = args->ROS10Factor;
    df_ptr->factor_actv = args->CROSActThreshold;
    df_ptr->cros = args->AllowCROS;

    // Calculate parameters
    calculate_BP(df_ptr, coef, &mainstruct, &sndstruct, &headstruct, &flankstruct, &backstruct);

    /*  ROSs DEBUG!   */
    if(args->verbose){
        std::cout << "*********** ROSs debug ************" << std::endl;
        std::cout <<  "-------Input Structure--------" << std::endl;
        std::cout <<  "fueltype: " << df_ptr->fueltype << std::endl;
        std::cout <<  "ws: " << df_ptr->ws << std::endl;
        std::cout <<  "waz: " << df_ptr->waz << std::endl;
        std::cout <<  "ps: " << df_ptr->ps << std::endl;
        std::cout <<  "saz: " << df_ptr->saz << std::endl;
        std::cout <<  "cur: " << df_ptr->cur << std::endl;
        std::cout <<  "elev: " << df_ptr->elev << std::endl;
        std::cout <<  "\n-------Mainout Structure--------" << std::endl;
        std::cout << "rss: " << mainstruct.rss << std::endl;
        std::cout << "angle: " << mainstruct.angle << std::endl;
        std::cout << "a:" << mainstruct.a << std::endl;
        std::cout << "b:" << mainstruct.b << std::endl;
        std::cout << "c:" << mainstruct.c << std::endl;
        std::cout << "covertype: " << mainstruct.covertype << std::endl;
        std::cout << "cros: " << mainstruct.cros << std::endl;
        std::cout <<  "\n-------Headout Structure--------" << std::endl;
        std::cout <<  "ros: " << headstruct.ros * args->HFactor << std::endl;
        std::cout <<  "rss: " << headstruct.rss << std::endl;
        std::cout <<  "\n------- Flank Structure--------" << std::endl;
        std::cout <<  "ros: " << flankstruct.ros * args->FFactor<< std::endl;
        std::cout <<  "rss: " << flankstruct.rss << std::endl;
        std::cout <<  "\n------- Back Structure--------" << std::endl;
        std::cout <<  "ros: " << backstruct.ros * args->BFactor << std::endl;
        std::cout <<  "rss: " << backstruct.rss << std::endl;
        std::cout <<  "\n------- Extra --------" << std::endl;
        std::cout <<  "lb: " << sndstruct.lb * args->BFactor << std::endl;
        std::cout << "*********** ROSs debug ************" << std::endl;
    }

//     double cartesianAngle = 270 - wdf_ptr->waz; 
//     if (cartesianAngle < 0){
//         cartesianAngle += 360;
//     } 

//     // Adjusting from Spanish forests angle
//     cartesianAngle =  wdf_ptr->waz;
//     double offset = cartesianAngle + 270;
//     cartesianAngle = 360 - (offset >= 360) * (cartesianAngle - 90) - (offset < 360) * offset;
//     if (cartesianAngle == 360)  cartesianAngle = 0;
//     if (cartesianAngle < 0) cartesianAngle += 360; 

    // Cartesian angle from Weather file
    double cartesianAngle = wdf_ptr->waz; 
    
    // ROS uncertainty
    double ROSRV = 0;
    if (args->ROSCV > 0) {
        ROSRV = randomROS;
    }

    // Display if verbose True (FBP ROSs, Main angle, and ROS std (if included))
    if (args->verbose) {
        std::cout << "Main Angle (raz): " << wdf_ptr->waz  << " Cartesian: " << cartesianAngle << std::endl;
        std::cout << "EllipseFactors : " << EllipseFactors[0] << ", " << EllipseFactors[1] << ", " << EllipseFactors[2] << std::endl;
        std::cout << "FBP Front ROS Value: " << headstruct.ros << std::endl; 
        std::cout << "FBP Flanks ROS Value: " << flankstruct.ros << std::endl;
        std::cout <<  "FBP Rear ROS Value: " << backstruct.ros << std::endl;
        std::cout << "FBP Front ROS Value: " << headstruct.ros * EllipseFactors[0] << std::endl; 
        std::cout << "FBP Flanks ROS Value: " << flankstruct.ros * EllipseFactors[1] << std::endl;
        std::cout <<  "FBP Rear ROS Value: " << backstruct.ros * EllipseFactors[2] << std::endl;
        std::cout <<  "EFactor Value: " << backstruct.ros * EllipseFactors[3] << std::endl;
        std::cout << "Std Normal RV for Stochastic ROS CV: " << ROSRV << std::endl;
    }

    // If cell cannot send (thresholds), then it will be burned out in the main loop
    double HROS = (1 + args->ROSCV * ROSRV) * headstruct.ros * EllipseFactors[0];

    // Extra debug step for sanity checks
    if (args->verbose){
            std::cout << "\nSending message conditions" << std::endl;
            std::cout << "HROS: " << HROS << " Threshold: "<<  args->ROSThreshold << std::endl;
    }

    // Check thresholds for sending messages
    if (HROS > args->ROSThreshold) {
        // True = -100
        repeat = -100;

        if (args->verbose) {
            std::cout << "\nRepeat condition: " << repeat << std::endl;
            std::cout << "Cell can send messages" << std::endl;
        }

        // ROS distribution method
        // Get wind
        double wss = wdf_ptr->ws;
        
        // Vector of k_opts
        // MHK TO DO: Add flag if KFactors 1-5 exist and create default values if not existing
        std::vector<double> k_opts;
        if (args->KFactor1) {
            k_opts.push_back(args->KFactor1);
            k_opts.push_back(args->KFactor2);
            k_opts.push_back(args->KFactor3);
            k_opts.push_back(args->KFactor4);
            k_opts.push_back(args->KFactor5);
        }

        // else {
        //     k_opts.push_back(1);
        //     k_opts.push_back(2);
        //     k_opts.push_back(0);
        //     k_opts.push_back(0);
        //     k_opts.push_back(0);            
        // }

        // ROS Boosters (placeholder in case they are needed in the future, for each neighbor)
        std::vector<double> ros_boosters;
        ros_boosters.push_back(1.);
        ros_boosters.push_back(1.);
        ros_boosters.push_back(1.);
        ros_boosters.push_back(1.);
        ros_boosters.push_back(1.);
        ros_boosters.push_back(1.);
        ros_boosters.push_back(1.);
        ros_boosters.push_back(1.);
        
        // Declare variable obtain dereferenced pointer (containing dictionary of K opt values)
        std::unordered_map<double, std::vector<double>> koptDict = *kdict_ptr;

        // Distribute ROS
        ros_distr(cartesianAngle,
                  headstruct.ros * EllipseFactors[0],
                  flankstruct.ros * EllipseFactors[1], 
                  backstruct.ros * EllipseFactors[2],
                  sndstruct.lb,
                  wss,
                  args->EllipticalOption,
                  args->CenteredDistribution,
                  args->NoAngleOffset,
                  EllipseFactors[3],
                  kdict_ptr, // From MangeFire --> ros_distr --> ellipse_adjusted (MHK)
                  k_opts,
                  ros_boosters,
                  args->LBFormula,
                  args,
                //   args->KOption,
                  args->verbose);

        // this is a iterator through the keyset of a dictionary
        for (auto&  _angle : this->ROSAngleDir) {
            double angle = _angle.first;
            int nb = angleToNb[angle];
            double ros = (1 + args->ROSCV * ROSRV) * _angle.second;

            if (args->verbose) {
                std::cout << "     (angle, realized ros in m/min): (" << angle << ", " << ros << ")" << std::endl;
            }

            // Workaround PeriodLen in 60 minutes
            this->fireProgress[nb] += ros * args->FirePeriodLen;   // Updates fire progress

            // If the message arrives to the adjacent cell's center, send a message
            if (this->fireProgress[nb] >= this->distToCenter[nb]) {
                msg_list.push_back(nb);
                FSCell->push_back(double(this->realId));
                FSCell->push_back(double(nb));
                FSCell->push_back(double(period));
                FSCell->push_back(ros);
                // FSCell->push_back(angle);
            }    

            // Info for debugging status of the cell and fire evolution			
            if (this->fireProgress[nb] < this->distToCenter[nb] && repeat == -100 && -100  != msg_list_aux[0]){
                    if (args->verbose){
                        std::cout << "A Repeat = TRUE flag is sent in order to continue with the current fire....." << std::endl;
                        std::cout << "Main workaround of the new sim logic....." << std::endl;
                    }
                    msg_list_aux[0] = repeat;
            }

        }

        if (args->verbose){
            printf("fireProgress Dict: ");
            for (auto & nb : this->fireProgress){
                std::cout << " " << nb.first << " : " << nb.second;
            }
            std::cout << std::endl;
        }
    }


    // If original is empty (no messages but fire is alive if aux_list is not empty)
    if  (msg_list.size() == 0){
        if (msg_list_aux[0] == -100){
            msg_list = msg_list_aux;
        }

        else{
            this->status = 2;   // we are done sending messages, call us burned
        }
    }

    if (args->verbose){
        std::cout << " ----------------- End of new manageFire function -----------------" << std::endl;
    }
    return msg_list;
}



/********************************************************************************************
    Get burned new logic: Checks if the ROS on its side is above a threshold for burning
        
    Returns     boolean  
    
    Inputs:
        period      int
        NMsg        int
        Season      int
        verbose     boolean
        df          Data frame
        coef        pointer
        args        hashmap    
********************************************************************************************/
bool CellsFBP::get_burned(int period, 
                          int season,
                          int NMsg,
                          inputs df[],  
                          fuel_coefs * coef,
                        //   koptDF * kopt_ptr,
                          arguments * args,
                          weatherDF * wdf_ptr) {
    if (args->verbose) { 
        std::cout << "ROS Threshold get_burned method" << std::endl;
        std::cout << "ROSThreshold: " << args->ROSThreshold << std::endl;
    }

    // Structures
    main_outs mainstruct;
    snd_outs sndstruct;
    fire_struc headstruct, backstruct, flankstruct;

    // Compute main angle and ROSs: forward, flanks and back
    df[this->id].waz = wdf_ptr->waz;
    df[this->id].ws = wdf_ptr->ws;
    df[this->id].mc1 = wdf_ptr->mc1;
    df[this->id].mc10 = wdf_ptr->mc10;
    df[this->id].mc100 = wdf_ptr->mc100;
    df[this->id].mcWoody = wdf_ptr->mcWoody;
    df[this->id].mcHerb = wdf_ptr->mcHerb;
    df[this->id].factor_cbd = args->CBDFactor;   
    df[this->id].factor_ccf = args->CCFFactor;
    df[this->id].factor_ros10 = args->ROS10Factor;
    df[this->id].factor_actv = args->CROSActThreshold;
    df[this->id].cros = args->AllowCROS;

    // Calculate
    // calculate(&(df[this->id]), coef, &mainstruct, &sndstruct, &headstruct, &flankstruct, &backstruct);
    calculate_BP(&(df[this->id]), coef, &mainstruct, &sndstruct, &headstruct, &flankstruct, &backstruct);

    if (args->verbose) { 
        std::cout << "\nMain Angle :" << wdf_ptr->waz << std::endl;
        std::cout << "Front ROS Value :" << headstruct.rss * args->HFactor << std::endl;
        std::cout << "Flanks ROS Value :" << flankstruct.rss * args->FFactor << std::endl;
        std::cout << "Rear ROS Value :" << backstruct.rss * args->BFactor << std::endl;
    }

    // Check a threshold for the ROS
    if (headstruct.rss  * args->HFactor > args->ROSThreshold) {
        this->status = 1;
        this->fireStarts = period;
        this->fireStartsSeason = season;
        this->burntP = period;
        return true;
    }
    // Not burned
    return false; 
}


/********************************************************************************************
                                    Old functions
********************************************************************************************/
/*
    Returns            void
    
    Inputs:  
    AdjacentCells      dictionary{string:[array integers]}
*/
void CellsFBP::set_Adj(std::unordered_map<std::string, int> & adjacentCells) {   
    // TODO: in python, these are pointers, maybe make these pointers too :P
    this->adjacents = adjacentCells;
}


/* 
    Returns            void
    
    Inputs:  
    Status_int         int
*/
void CellsFBP::setStatus(int status_int) {  
    this->status = status_int;
}


/*
    Returns            string
    
    Inputs:  
*/
std::string CellsFBP::getStatus() {
    // Return cell's status
    return this->StatusD[this->status];
}


/*
    Returns           boolean
    
    Inputs:
    period            int
    Season            int
    IgnitionPoints    array of int
    df                Data frame
    coef              pointer
    ROSThresh         double
    HFIThreshold      double
 */
bool CellsFBP::ignition(int period, 
                        int year,
                        std::vector<int> & ignitionPoints,
                        inputs * df_ptr,
                        fuel_coefs * coef,
                        // koptDF * kopt_ptr,
                        arguments *args,
                        weatherDF * wdf_ptr) {
    
    // If we have ignition points, update
    if (ignitionPoints.size() > 0) {
        this->status = 1;
        this->fireStarts = period;
        this->fireStartsSeason = year;
        this->burntP = period;

        // An ignition has happened
        return true;
    } else {
        // Ignites if implied head ros andfire intensity are high enough
        main_outs mainstruct;
        snd_outs sndstruct;
        fire_struc headstruct, backstruct, flankstruct;

        // Populate inputs 
        df_ptr->waz = wdf_ptr->waz;
        df_ptr->ws = wdf_ptr->ws;
        // df_ptr->scen = wdf_ptr->scenario;
        df_ptr->factor_cbd = args->CBDFactor;   
        df_ptr->factor_ccf = args->CCFFactor;
        df_ptr->factor_ros10 = args->ROS10Factor;
        df_ptr->factor_actv = args->CROSActThreshold;
        df_ptr->cros = args->AllowCROS;

        // calculate(df_ptr, coef, &mainstruct, &sndstruct, &headstruct, &flankstruct, &backstruct);
        calculate_BP(df_ptr, coef, &mainstruct, &sndstruct, &headstruct, &flankstruct, &backstruct);

        if (args->verbose) {
            std::cout << "\nIn ignition function" << std::endl;
            std::cout << "Main Angle: " << wdf_ptr->waz << std::endl;
            std::cout << "Front ROS Value: " << headstruct.rss * args->HFactor << std::endl;
            std::cout << "Flanks ROS Value: " << flankstruct.rss * args->FFactor << std::endl;
            std::cout << "Rear ROS Value: " << backstruct.rss * args->BFactor << std::endl;
        }

        // Check a threshold for the ROS
        if (headstruct.rss * args->HFactor > args->ROSThreshold) {
            if (args->verbose) {
                std::cout << "Head (ROS, FI) values of: (" << headstruct.rss * args->HFactor<< ", " << ") are enough for ignition" << std::endl;
            }

            this->status = 1;
            this->fireStarts = period;
            this->fireStartsSeason = year;
            this->burntP = period;

            return true;
        }
        return false;
    }
}



/*
    Returns      void
    Inputs
    ID           int
    period       int
*/
void CellsFBP::harvested(int id, int period) {     
    // TODO: unused param
    this->status = 3;
    this->harvestStarts = period;
}



/*
    Returns      void
*/
void CellsFBP::print_info() {    
    std::cout << "Cell Information" << std::endl;
    std::cout << "ID = "  << this->id<< std::endl;
    std::cout << "In Forest ID = "  << this->realId<< std::endl;
    std::cout << "Status = " << this->StatusD[this->status] << std::endl;
    std::cout << "Coordinates: ";
    std::cout << this->coord[0] <<  " " << this->coord[1]  << std::endl;

    std::cout << "Area = "<<  this->area << std::endl;
    std::cout << "FTypes = "<< this->FTypeD[this->fType] << std::endl;
    std::cout << "AdjacentCells:";
    for (auto & nb : this->adjacents){
        std::cout << " " << nb.first << ":" << nb.second;
    }
    std::cout << std::endl;

    printf("Angle Dict: ");
    for (auto & nb : this->angleDict){
        std::cout << " " << nb.first << " : " << nb.second;
    }
    std::cout << std::endl;

    printf("Ros Angle Dict: ");
    for (auto & nb : this->ROSAngleDir){
        std::cout << " " << nb.first << " : " << nb.second;
    }
    std::cout << std::endl;


    printf("angleToNb Dict: ");
    for (auto & nb : this->angleToNb){
        std::cout << " " << nb.first << " : " << nb.second;
    }
    std::cout << std::endl;


    printf("fireProgress Dict: ");
    for (auto & nb : this->fireProgress){
        std::cout << " " << nb.first << " : " << nb.second;
    }
    std::cout << std::endl;


    printf("distToCenter Dict: ");
    for (auto & nb : this->distToCenter){
        std::cout << " " << nb.first << " : " << nb.second;
    }
    std::cout << std::endl;
}