// Definition
#ifndef CELLSFBP
#define CELLSFBP

// Basic libraries
#include "fuelmodelBP.h"
#include "ReadCSV.h"
#include "ReadArgs.h"
#include "Ellipse.h"
#include <stdio.h>
#include <math.h>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <vector>

// Torch
// #include <torch/torch.h>
// #include <torch/script.h> 

using namespace std;

class CellsFBP {
    // TODO: find where to put the enums
    public:
        // immutable [cells]
        int id;
        int fType;
        int realId;
        double _ctr2ctrdist;
        double area;
        double perimeter;
        
        // Conditions
        int nftype;
        float elev ;
        float saz ;
        float ws ;
        float waz ;
        int fm_1h ;
        int fm_10h ;
        int fm_100h ;
        int fm_dead ;
        int fm_live ;

        // Extras
        std::string fType2;
        std::vector<int> coord; //maybe change this into a tuple or class	CP: 2-tuple (int)
        std::unordered_map<std::string, int> adjacents; // CP: dictionary {string: [int array]}
        
        // Read and load k_opt table from csv
        std::unordered_map<double, std::vector<double>> koptDict; // keys = double, values = double
        std::unordered_map<double, std::vector<double>> *kdict_ptr = new std::unordered_map<double, std::vector<double>>;


        // Ftype dicts
        string FTypeD[3];
        string StatusD[5];

        // mutable
        int status;
        int hPeriod;
        int fireStarts;
        int harvestStarts;
        int fireStartsSeason;
        int burntP;
        int tYears;	

        std::unordered_map<int, std::vector<int>> gMsgList; // {40 -> [1, 2, 3] }
        std::unordered_map<int, std::vector<int>> gMsgListSeason;
        std::unordered_map<int, double> fireProgress;      // CP: dictionary {int: double}
        std::unordered_map<int, double> angleDict;         // CP: dictionary {int: double}
        std::unordered_map<int, double> ROSAngleDir;       // CP: dictionary {int: double|None}  
        std::unordered_map<int, double> distToCenter;      // CP: dictionary {int: double}
        std::unordered_map<int, int> angleToNb;            // CP: dictionary {double: int}

        // TODO: reference to shared object

        // Constructor and methods here
        CellsFBP(int _id, 
                 double _area,
                 std::vector<int> _coord, 
                 int _fType, 
                 std::string _fType2, 
                 double _perimeter, 
                 int _status, std::unordered_map<std::string,
                 int> & _adjacents, 
                 int _realId);

        // Fire fields init
        void initializeFireFields(std::vector<std::vector<int>> & coordCells,
                                  std::unordered_set<int> & availSet); 

        // Elliptical fits
        std::vector<double> alexander_ellipse(double ros_mmin, 
                                              double ws_kmhr,
                                              int lb_estimation_scheme,
                                              bool verbose);
        std::vector<double> fbp_ellipse(double forward,
                                        double flank,
                                        double back,
                                        bool verbose);
        std::vector<double> lb_ellipse(double ros_mmin, 
                                       double LB,
                                       bool verbose);
        std::vector<double> ellipse_adjusted(double hros_mmin,
                                             double ws_kmhr,
                                             std::vector<double> k_opts,
                                             std::unordered_map<double, std::vector<double>> *kdict_ptr,
                                             int lb_estimation_scheme,
                                             std::string KOption,
                                             bool verbose);
    
        // ROS calculations
        void ros_distr_old(double thetafire, 
                           double forward,
                           double flank,
                           double back);

        double rhoTheta(double theta, 
                        double a,
                        double b);
        double rhoTheta_center(double theta, 
                               double a,
                               double b);
        double rhoTheta_adjusted(double angle,
                                 double LW,
                                 double Ebar,
                                 double Ebar_c,
                                 double hros_mmin_c,
                                 std::vector<double> k_opts,
                                 bool deactivate_angle_offset,
                                 bool verbose);
    
    
        // ROS distribution
        void ros_distr(double thetafire,
                       double forward,
                       double flank,
                       double back,
                       double lb,
                       double ws,
                       int elliptical_option,
                       bool centered_distribution,
                       bool deactivate_angle_offset,
                       double EFactor,
                       std::unordered_map<double, std::vector<double>> *kdict_ptr,
                       std::vector<double> k_opts,
                       std::unordered_map<int, double> ros_boosters,
                       int lb_estimation_scheme,
                       arguments * args,
                       bool verbose);
    
        // General LB
        double general_lb(double ws_kmhr,
                          int lb_estimation_scheme);
    
    
        // Manage Fire
        std::vector<int> manageFire(int period,
                                    std::unordered_set<int> & AvailSet,      
                                    inputs df[], fuel_coefs * coef, 
                                    std::unordered_map<double, std::vector<double>> *kdict_ptr,
                                    std::vector<std::vector<int>> & coordCells, 
                                    std::unordered_map<int, CellsFBP> & Cells_Obj, 
                                    arguments * args, 
                                    weatherDF * wdf_ptr,
                                    std::vector<double> * FSCell,
                                    double randomROS);

        std::vector<int> manageFireBBO(int period,
                                       std::unordered_set<int> & AvailSet,      
                                       inputs * df_ptr, fuel_coefs * coef, 
                                       std::unordered_map<double, std::vector<double>> *kdict_ptr,
                                       std::vector<std::vector<int>> & coordCells, 
                                       std::unordered_map<int, CellsFBP> & Cells_Obj, 
                                       arguments * args,
                                       weatherDF * wdf_ptr,
                                       std::vector<double> * FSCell,
                                       double randomROS, 
                                       std::vector<float> & EllipseFactors);

        // Get burned conditions
        bool get_burned(int period, 
                        int season,
                        int NMsg,
                        inputs df[],
                        fuel_coefs * coef, 
                        // koptDF * kopt_ptr,
                        // koptDict * kdict_ptr,
                        arguments * args,
                        weatherDF * wdf_ptr);

        // Extras
        void set_Adj(std::unordered_map<std::string, int> & adjacentCells);
        void setStatus(int status_int);
        std::string getStatus();
        bool ignition(int period,
                      int year,
                      std::vector<int> & ignitionPoints,
                      inputs * df_ptr, 
                      fuel_coefs * coef,
                    //   koptDF * kopt_ptr,
                    //   koptDict * kdict_ptr,
                      arguments *args,
                      weatherDF * wdf_ptr);
        void harvested(int id, int period);
        void print_info();
        double round_up(double value, 
                        int decimal_places);

    
    // Private methods
    private:
        double allocate(double offset,
                        double base,
                        double ros1,
                        double ros2);
        double slope_effect(double elev_i,
                            double elev_j,
                            double cellsize);
};

#endif
