#ifndef CELLS
#define CELLS

// include stuff
#include "Ellipse.h"
#include "ReadArgs.h"

#include <math.h>
#include <stdio.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;

std::vector<int> adjacentCells(int cell, int nrows, int ncols);
/*
 *   Weather structure
 */
typedef struct
{
    float ws, waz, rh, tmp, apcp, ffmc, dmc, dc, isi, bui,
        fwi;  // David: some variables only used on C2FK and not on C2FSB and
              // viceversa
} weatherDF;  // David: Moved here to simplify inclusion

typedef struct
{
    char fueltype[4];
    float ws, saz, cur, ps, cbh, ccf, cbd, elev, tmp, rh, lat, lon, ffmc, bui, gfl,
        tree_height;  // David: some variables only used on C2FK and not on C2FSB and
              // viceversa
    int waz, nftype, FMC, time, pattern, mon, jd, jd_min, pc, pdf;
} inputs;  // David: Moved here to simplify inclusion

typedef struct
{
    char fueltype[4];
    float p1, p2, p3;    // hros coef
    float q1, q2, q3;    // flame length coef
    float q, bui0, cfl;  // fbp params
    float cbh, fmc, fl, h;
    int nftype;
} fuel_coefs;

typedef struct
{
    float hffmc, sfc, csi, fl, fh, a, b, c, rss, angle, ros_active, cfb, se, rso, fmc, sfi, isi, be, sf, raz, wsv, ff,
        crown_intensity, crown_flame_length, max_flame_length;
    char covertype;
    int crown, jd_min, jd;
} main_outs;

typedef struct
{
    float ros, dist, rost, cfb, fc, cfc, time, rss, isi;
    char fd;
    double fi;
} fire_struc;

typedef struct
{
    float lb, area, perm, pgr, lbt;
} snd_outs;

class Cells
{
    // TODO: find where to put the enums
  public:
    // immutable
    int id;
    int fType;
    int realId;
    double _ctr2ctrdist;
    double area;
    double perimeter;

    std::string fType2;
    std::vector<int> coord;  // maybe change this into a tuple or class	CP: 2-tuple (int)
    // std::unordered_map<std::string, int> adjacents; // CP: dictionary {string:
    // [int array]}

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

    std::unordered_map<int, std::vector<int>> gMsgList;  // {40 -> [1, 2, 3] }
    std::unordered_map<int, std::vector<int>> gMsgListSeason;
    std::unordered_map<int, double> fireProgress;  // CP: dictionary {int: double}
    std::unordered_map<int, double> angleDict;     // CP: dictionary {int: double}
    std::unordered_map<int, double> ROSAngleDir;   // CP: dictionary {int: double|None}   Instead of None we
                                                   // can use a determined number like -9999 = None  TODO:
                                                   // maybe int : double
    std::unordered_map<int, double> distToCenter;  // CP: dictionary {int: double}
    std::unordered_map<int, int> angleToNb;        // CP: dictionary {double: int}

    // TODO: reference to shared object

    // constructor and methods here
    Cells(int _id,
          double _area,
          std::vector<int> _coord,
          int _fType,
          std::string _fType2,
          double _perimeter,
          int _status,
          int _realId);

    void initializeFireFields(std::vector<std::vector<int>>& coordCells,
                              std::unordered_set<int>& availSet,
                              int cols,
                              int rows);  // TODO: need TYPE

    // Elliptical fits (v2)
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
    
    // General LB (v2)
    double general_lb(double ws_kmhr,
                        int lb_estimation_scheme);

    // Read and load k_opt table from csv (K-factor optimization)
    std::unordered_map<double, std::vector<double>> koptDict; // keys = double, values = double
    std::unordered_map<double, std::vector<double>> *kdict_ptr = new std::unordered_map<double, std::vector<double>>;

    void ros_distr_old(double thetafire, double forward, double flank, double back);
    double rhoTheta(double theta, double a, double b);
    void ros_distr(double thetafire, double forward, double flank, double back, double EFactor);
    void ros_distr_V2(double thetafire, double a, double b, double c, double EFactor);
    
    // ROS distribution for US (v2)
    void ros_distr_US(double thetafire,
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
    double rhoTheta_adjusted(double angle,
                            double LW,
                            double Ebar,
                            double Ebar_c,
                            double hros_mmin_c,
                            std::vector<double> k_opts,
                            bool deactivate_angle_offset,
                            bool verbose);
    double rhoTheta_center(double theta, 
                            double a,
                            double b);                            

    std::vector<int> manageFire(int period,
                                std::unordered_set<int>& AvailSet,
                                inputs df[],
                                fuel_coefs* coef,
                                std::vector<std::vector<int>>& coordCells,
                                std::unordered_map<int, Cells>& Cells_Obj,
                                arguments* args,
                                weatherDF* wdf_ptr,
                                std::vector<double>* FSCell,
                                std::vector<float>* crownMetrics,
                                bool& activeCrown,
                                double randomROS,
                                int perimeterCells,
                                std::vector<int>& crownState,
                                std::vector<float>& crownFraction,
                                std::vector<float>& surfFraction,
                                std::vector<float>& Intensities,
                                std::vector<float>& RateOfSpreads,
                                std::vector<float>& SurfaceFlameLengths,
                                std::vector<float>& CrownFlameLengths,
                                std::vector<float>& CrownIntensities,
                                std::vector<float>& MaxFlameLengths);

    std::vector<int> manageFireBBO(int period,
                                   std::unordered_set<int>& AvailSet,
                                   inputs* df_ptr,
                                   fuel_coefs* coef,
                                   std::vector<std::vector<int>>& coordCells,
                                   std::unordered_map<int, Cells>& Cells_Obj,
                                   arguments* args,
                                   weatherDF* wdf_ptr,
                                   std::vector<double>* FSCell,
                                   std::vector<float>* crownMetrics,
                                   bool& activeCrown,
                                   double randomROS,
                                   int perimeterCells,
                                   std::vector<float>& EllipseFactors,
                                   std::vector<int>& crownState,
                                   std::vector<float>& crownFraction,
                                   std::vector<float>& surfFraction,
                                   std::vector<float>& Intensities,
                                   std::vector<float>& RateOfSpreads,
                                   std::vector<float>& FlameLengths);

    bool get_burned(int period,
                    int season,
                    int NMsg,
                    inputs df[],
                    fuel_coefs* coef,
                    arguments* args,
                    weatherDF* wdf_ptr,
                    bool& activeCrown,
                    int perimeterCells);

    // void set_Adj(std::unordered_map<std::string, int> & adjacentCells);

    void setStatus(int status_int);

    std::string getStatus();

    bool ignition(int period,
                  int year,
                  std::vector<int>& ignitionPoints,
                  inputs* df_ptr,  // WORKING CHECK OK
                  fuel_coefs* coef,
                  arguments* args,
                  weatherDF* wdf_ptr,
                  bool& activeCrown,
                  int perimeterCells);

    void harvested(int id, int period);

    void print_info();

    // Utility (v2)
    double round_up(double value, int decimal_places); 

  private:
    double allocate(double offset, double base, double ros1, double ros2);
    float slope_effect(float elev_i, float elev_j, int cellsize);
};

#endif
