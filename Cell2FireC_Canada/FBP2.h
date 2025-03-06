#ifndef FBP2
#define FBP2

// Define your structs
typedef struct {
    char fueltype[4];
    float ffmc, ws, gfl, bui, lat, lon;
    int time, pattern, mon, jd, jd_min, waz, ps, saz, pc, pdf, cur, elev, hour, hourly;
} inputs;

typedef struct {
    char fueltype[4];
    float q, bui0, cbh, cfl;
    double a, b, c;
} fuel_coefs;

typedef struct {
    float ros, dist, rost, cfb, fc, cfc, time, rss, isi;
    char fd;
    double fi;
} fire_struc;

typedef struct {
    float hffmc, sfc, csi, rso, fmc, sfi, rss, isi, be, sf, raz, wsv, ff;
    int jd_min, jd;
    char covertype;
} main_outs;

typedef struct {
    float lb, area, perm, pgr, lbt;
} snd_outs;

// Declare your functions
float rate_of_spread(inputs *inp, fuel_coefs *ptr, main_outs *at);
float ros_calc(inputs *inp, fuel_coefs *ptr, float isi, float *mu);
float bui_effect(fuel_coefs *ptr, main_outs *at, float bui);
float ISF_mixedwood(fuel_coefs *ptr, float isz, int pc, float sf);
float ISF_deadfir(fuel_coefs *ptr, float isz, int pdf, float sf);
float grass(fuel_coefs *ptr, float cur, float isi, float *mu);
float mixed_wood(fuel_coefs *ptr, float isi, float *mu, int pc);
float dead_fir(fuel_coefs *ptr, int pdf, float isi, float *mu);
float D2_ROS(fuel_coefs *ptr, float isi, float bui, float *mu);
float slope_effect(inputs *inp, fuel_coefs *ptr, main_outs *at, float isi);
float backfire_isi(main_outs *at);
float flankfire_ros(float ros, float bros, float lb);

#endif // YOUR_HEADER_FILE_H
