#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstring>

// Structure definitions
struct inputs {
    float ps;
    float ffmc;
    float ws;
    float waz;
    float saz;
    float cur;
    float bui;
    int pc;
    int pdf;
    char fueltype[4];
};

struct fuel_coefs {
    float a;
    float b;
    float c;
    float q;
    float bui0;
    float cbh;
    float cfl;
    char fueltype[4];
};

struct main_outs {
    float ff;
    float wsv;
    float isi;
    float rss;
    float sf;
    float raz;
    float be;
};

// Function prototypes
float rate_of_spread(inputs *inp, fuel_coefs *ptr, main_outs *at);
float ffmc_effect(float ffmc);
float slope_effect(inputs *inp, fuel_coefs *ptr, main_outs *at, float isi);
float ISF_mixedwood(fuel_coefs *ptr, float isz, int pc, float sf);
float ISF_deadfir(fuel_coefs *ptr, float isz, int pdf, float sf);
float ros_calc(inputs *inp, fuel_coefs *ptr, float isi, float *mult);
float conifer(fuel_coefs *ptr, float isi, float *mu);
float grass(fuel_coefs *ptr, float cur, float isi, float *mult);
float mixed_wood(fuel_coefs *ptr, float isi, float *mu, int pc);
float dead_fir(fuel_coefs *ptr, int pdf, float isi, float *mu);
float D2_ROS(fuel_coefs *ptr, float isi, float bui, float *mu);
float bui_effect(fuel_coefs *ptr, main_outs *at, float bui);
float backfire_ros(inputs *inp, fuel_coefs *ptr, main_outs *at, float bisi);
float flankfire_ros(float ros, float bros, float lb);

// Global variables
float slopelimit_isi = 0.01;
float numfuels = 18;

// Function implementations
float rate_of_spread(inputs *inp, fuel_coefs *ptr, main_outs *at)
{
float fw,isz,mult,*mu=&mult,rsi  ;
at->ff=ffmc_effect(inp->ffmc);
at->raz=inp->waz;
isz=0.208*at->ff;
if(inp->ps>0)at->wsv=slope_effect(inp,ptr,at,isz);
else at->wsv=inp->ws;
if(at->wsv<40.0)fw=exp(0.05039*at->wsv);
else fw=12.0*(1.0-exp(-0.0818*(at->wsv-28)));
at->isi=isz*fw;
rsi=ros_calc(inp,ptr,at->isi,mu);
at->rss=rsi*bui_effect(ptr,at,inp->bui);
return(at->rss);
}

float ffmc_effect(float ffmc)
{
float mc,ff;
mc=147.2*(101.0-ffmc)/(59.5+ffmc);
ff=91.9*exp(-0.1386*mc)*(1+pow(mc,5.31)/49300000.0);
return (ff);
}

float slope_effect(inputs *inp,fuel_coefs *ptr,main_outs *at, float isi )
/* ISI is ISZ really */
{
double isf,rsf,wse,ps,rsz,wsx,wsy,wsex,wsey,wsvx,wsvy,
    wrad,srad,wsv,raz,check,wse2,wse1;
float mu=0.0;
ps=inp->ps*1.0;

if(ps>70.0)ps=70.0;   /* edited in version 4.6*/
at->sf=exp(3.533*pow(ps/100.0,1.2));
if(at->sf>10.0)at->sf=10.00;  /* added to ensure maximum is correct in version 4.6  */

if(strncmp(ptr->fueltype,"M1",2)==0 || strncmp(ptr->fueltype,"M2",2)==0)
    isf=ISF_mixedwood(ptr,isi,inp->pc,at->sf);
else if(strncmp(ptr->fueltype,"M3",2)==0 || strncmp(ptr->fueltype,"M4",2)==0)
    isf=ISF_deadfir(ptr,isi,inp->pdf,at->sf);
else{
    rsz=ros_calc(inp,ptr,isi,&mu);
    rsf=rsz*at->sf;

    if(rsf>0.0)check=1.0-pow((rsf/(mu*ptr->a)),(1.0/ptr->c) );
    else check=1.0;
    if(check<slopelimit_isi)check=slopelimit_isi;

    isf=(1.0/(-1.0*ptr->b))*log(check);
}
if(isf==0.0)isf=isi;  /* should this be 0.0001 really  */
wse1 = log(isf/(0.208*at->ff))/0.05039;
if(wse1<=40.0) wse=wse1;
else{
    if(isf>(0.999*2.496*at->ff) ) isf=0.999*2.496*at->ff;
    wse2=28.0-log(1.0-isf/(2.496*at->ff))/0.0818;
    wse=wse2;
}
wrad=inp->waz/180.0*3.1415926;
wsx=inp->ws*sin(wrad);
wsy=inp->ws*cos(wrad);
srad=inp->saz/180.0*3.1415926;
wsex=wse*sin(srad);
wsey=wse*cos(srad);
wsvx=wsx+wsex;
wsvy=wsy+wsey;
wsv=sqrt(wsvx*wsvx+wsvy*wsvy);
raz=acos(wsvy/wsv);
raz=raz/3.1415926*180.0;
if(wsvx<0)raz=360-raz;
at->raz=raz;
return( (float)(wsv) );
}

float ISF_mixedwood(fuel_coefs *ptr, float isz, int pc, float sf)
{
float check, mult,rsf_d1,rsf_c2,isf_d1,isf_c2;
int i;

rsf_c2=sf*ptr->a*pow( (1.0-exp(-1.0*ptr->b*isz)), ptr->c);
if(rsf_c2>0.0)check=1.0-pow((rsf_c2/(ptr->a)),(1.0/ptr->c) );
else check=1.0;
if(check<slopelimit_isi)check=slopelimit_isi;
isf_c2=(1.0/(-1.0*ptr->b))*log(check);

if(strncmp(ptr->fueltype,"M2",2)==0) mult=0.2;
else mult=1.0;
for(i=0;strncmp(ptr->fueltype,"D1",2)!=0 && i<numfuels;ptr++,i++);
rsf_d1=sf*(mult*ptr->a)*pow( (1.0-exp(-1.0*ptr->b*isz) ),ptr->c);

if(rsf_d1>0.0)check=1.0-pow((rsf_d1/(mult*ptr->a)),(1.0/ptr->c) );
else check=1.0;
if(check<slopelimit_isi)check=slopelimit_isi;
isf_d1=(1.0/(-1.0*ptr->b))*log(check);

return  ( ((float)(pc)/100.0)*isf_c2 + (100-(float)(pc))/100.0*isf_d1  );
}

float ISF_deadfir(fuel_coefs *ptr, float isz, int pdf, float sf)
{
float check, mult,rsf_d1,rsf_max,isf_d1,isf_max;
int i;

rsf_max=sf*ptr->a*pow( (1.0-exp(-1.0*ptr->b*isz)), ptr->c);
if(rsf_max>0.0)check=1.0-pow((rsf_max/(ptr->a)),(1.0/ptr->c) );
else check=1.0;
if(check<slopelimit_isi)check=slopelimit_isi;
isf_max=(1.0/(-1.0*ptr->b))*log(check);

if(strncmp(ptr->fueltype,"M4",2)==0) mult=0.2;
else mult=1.0;

for(i=0;strncmp(ptr->fueltype,"D1",2)!=0 && i<numfuels;ptr++,i++);
rsf_d1=sf*(mult*ptr->a)*pow( (1.0-exp(-1.0*ptr->b*isz) ),ptr->c);

if(rsf_d1>0.0)check=1.0-pow((rsf_d1/(mult*ptr->a)),(1.0/ptr->c) );
else check=1.0;
if(check<slopelimit_isi)check=slopelimit_isi;
isf_d1=(1.0/(-1.0*ptr->b))*log(check);

return  ( ((float)(pdf)/100.0)*isf_max + ( 100.0-(float)(pdf) )/100.0*isf_d1  );
}

float conifer(fuel_coefs *ptr, float isi, float *mu)
{
*mu=1.0;
return( ptr->a*pow( (1.0-exp(-1.0*ptr->b*isi) ),ptr->c) );
}


float ros_calc(inputs *inp, fuel_coefs *ptr,float isi,float *mult)
{
float ros;
if(strncmp(inp->fueltype,"O1",2)==0)
            return ( grass(ptr, inp->cur ,isi,mult) );
if(strncmp(inp->fueltype,"M1",2)==0 || strncmp(inp->fueltype,"M2",2)==0)
            return ( mixed_wood(ptr,isi,mult,inp->pc) );
if(strncmp(inp->fueltype,"M3",2)==0 || strncmp(inp->fueltype,"M4",2)==0)
            return ( dead_fir(ptr,inp->pdf,isi,mult) );
if(strncmp(inp->fueltype,"D2",2)==0)
            return ( D2_ROS(ptr,isi,inp->bui,mult) );
/* if all else has fail its a conifer   */
return ( conifer(ptr,isi,mult));
}

float grass(fuel_coefs *ptr,float cur,float isi,float *mult)
{
float mu,ros;
if((float)(cur)>=58.8) mu=0.176 + 0.02*( (float)(cur) - 58.8 ) ;
else mu=0.005*(exp(0.061*(float)(cur) ) - 1.0) ;
ros=mu * (ptr->a*pow( (1.0-exp(-1.0*ptr->b*isi)) , ptr->c));
if(mu<0.001)mu=0.001;  /* to have some value here*/
*mult=mu;
return(ros);
}

float mixed_wood(fuel_coefs *ptr, float isi,float *mu, int pc)
{
float ros, mult,ros_d1,ros_c2;
int i;
*mu=pc/100.0;
ros_c2=ptr->a*pow( (1.0-exp(-1.0*ptr->b*isi)), ptr->c);
if(strncmp(ptr->fueltype,"M2",2)==0) mult=0.2;
else mult=1.0;
for(i=0;strncmp(ptr->fueltype,"D1",2)!=0 && i<numfuels;ptr++,i++);
if(i>=numfuels) { printf(" prob in mixedwood   d1 not found \n");exit(9);}
ros_d1=ptr->a*pow( (1.0-exp(-1.0*ptr->b*isi) ),ptr->c);

ros=(pc/100.0)*ros_c2 + mult* (100-pc)/100.0*ros_d1;
return(ros);
}

float dead_fir(fuel_coefs *ptr, int pdf, float isi, float *mu)
{
double a,b,c;
int i;
float ros,rosm3or4_max,ros_d1, greenness=1.0;
if(strncmp(ptr->fueltype,"M4",2)==0)greenness=0.2;

rosm3or4_max=ptr->a*pow( ( 1.0-exp(-1.0*ptr->b*isi)),ptr->c);

for(i=0;strncmp(ptr->fueltype,"D1",2)!=0 && i<numfuels;ptr++,i++);
if(i>=numfuels) { printf(" prob in mixedwood   d1 not found \n");exit(9);}
ros_d1=ptr->a*pow( (1.0-exp(-1.0*ptr->b*isi) ),ptr->c);

ros=(float)(pdf)/100.0*rosm3or4_max + (100.0-(float)(pdf))/100.0*greenness*ros_d1;

*mu=(float)(pdf)/100.0;

return(ros);
}

float D2_ROS(fuel_coefs *ptr, float isi, float bui, float *mu)
{
*mu=1.0;
if(bui>=80) return( ptr->a*pow( (1.0-exp(-1.0*ptr->b*isi) ),ptr->c) );
else return (0.0);
}

float bui_effect(fuel_coefs *ptr,main_outs *at, float bui)
{
float  bui_avg=50.0;

if(bui==0) bui=1.0;
at->be=exp(bui_avg*log(ptr->q)*( (1.0/bui)-(1.0/ptr->bui0) ) );
return (at->be);
}

float backfire_ros(inputs *inp,fuel_coefs *ptr,main_outs *at,float bisi)
{
float mult=0.0,bros;
bros=ros_calc(inp,ptr,bisi,&mult);
bros *= bui_effect(ptr,at,inp->bui);
return(bros);
}

float flankfire_ros(float ros,float bros,float lb)
{
    return  ( (ros+bros)/(lb*2.0) );
}

// Utility function to trim leading and trailing spaces
std::string trim(const std::string& str) {
    auto start = str.find_first_not_of(" ");
    auto end = str.find_last_not_of(" ");
    return (start == std::string::npos || end == std::string::npos) ? "" : str.substr(start, end - start + 1);
}

// Utility function to safely convert string to float
float safe_stof(const std::string& str) {
    try {
        return str.empty() ? 0.0f : std::stof(str);
    } catch (const std::invalid_argument&) {
        std::cerr << "Invalid argument for float conversion: " << str << std::endl;
        return 0.0f;
    }
}

int safe_stoi(const std::string& str) {
    try {
        return str.empty() ? 0 : std::stoi(str);
    } catch (const std::invalid_argument&) {
        std::cerr << "Invalid argument for int conversion: " << str << std::endl;
        return 0;
    }
}

// Function to parse CSV line and convert values
void parse_CSV(const std::string& line) {
    std::stringstream ss(line);
    std::string value;

    std::vector<float> values;
    while (std::getline(ss, value, ',')) {
        try {
            float number = safe_stof(value);
            values.push_back(number);
        } catch (const std::invalid_argument& e) {
            std::cerr << "Conversion failed: " << e.what() << std::endl;
            values.push_back(0.0f);  // Or handle it as appropriate
        } catch (const std::out_of_range& e) {
            std::cerr << "Number out of range: " << e.what() << std::endl;
            values.push_back(0.0f);  // Or handle it as appropriate
        }
    }

    // Process the values as needed
for (std::vector<float>::iterator it = values.begin(); it != values.end(); ++it) {
    float val = *it;
        std::cout << val << " ";
    }
    std::cout << std::endl;
}



// Function to get fuel coefficients (dummy implementation)
fuel_coefs get_fuel_coefs(const std::string &fueltype) {
    fuel_coefs coefs;
    strcpy(coefs.fueltype, fueltype.c_str());
    coefs.a = 0.1;
    coefs.b = 0.01;
    coefs.c = 0.5;
    coefs.q = 0.5;
    coefs.bui0 = 0.0;
    coefs.cbh = 0.0;
    coefs.cfl = 0.0;
    return coefs;
}

// Function to set up fuel coefficients
void setup_const(fuel_coefs *ptr) {
    // Fuel types array initialization
    // Fuel type 0
    strncpy(ptr->fueltype, "M1", 3); ptr->a = 110.0; ptr->b = 0.0282; ptr->c = 1.5;
    ptr->q = 0.80; ptr->bui0 = 50; ptr->cbh = 6; ptr->cfl = 0.80;
    // Fuel type 1
    ptr++; strncpy(ptr->fueltype, "M2", 3); ptr->a = 110.0; ptr->b = 0.0282; ptr->c = 1.5;
    ptr->q = 0.80; ptr->bui0 = 50; ptr->cbh = 6; ptr->cfl = 0.80;
    // Fuel type 2
    ptr++; strncpy(ptr->fueltype, "M3", 3); ptr->a = 120.0; ptr->b = 0.0572; ptr->c = 1.4;
    ptr->q = 0.80; ptr->bui0 = 50; ptr->cbh = 6; ptr->cfl = 0.80;
    // Fuel type 3
    ptr++; strncpy(ptr->fueltype, "M4", 3); ptr->a = 100.0; ptr->b = 0.0404; ptr->c = 1.48;
    ptr->q = 0.80; ptr->bui0 = 50; ptr->cbh = 6; ptr->cfl = 0.80;
    // Fuel type 4
    ptr++; strncpy(ptr->fueltype, "C1", 3); ptr->a = 90.0; ptr->b = 0.0649; ptr->c = 4.5;
    ptr->q = 0.90; ptr->bui0 = 72; ptr->cbh = 2; ptr->cfl = 0.75;
    // Fuel type 5
    ptr++; strncpy(ptr->fueltype, "C2", 3); ptr->a = 110.0; ptr->b = 0.0282; ptr->c = 1.5;
    ptr->q = 0.70; ptr->bui0 = 64; ptr->cbh = 3; ptr->cfl = 0.80;
    // Fuel type 6
    ptr++; strncpy(ptr->fueltype, "C3", 3); ptr->a = 110.0; ptr->b = 0.0444; ptr->c = 3.0;
    ptr->q = 0.75; ptr->bui0 = 62; ptr->cbh = 8; ptr->cfl = 1.15;
    // Fuel type 7
    ptr++; strncpy(ptr->fueltype, "C4", 3); ptr->a = 110.0; ptr->b = 0.0293; ptr->c = 1.5;
    ptr->q = 0.80; ptr->bui0 = 66; ptr->cbh = 4; ptr->cfl = 1.20;
    // Fuel type 8
    ptr++; strncpy(ptr->fueltype, "C5", 3); ptr->a = 30.0; ptr->b = 0.0697; ptr->c = 4.0;
    ptr->q = 0.80; ptr->bui0 = 56; ptr->cbh = 18; ptr->cfl = 1.20;
    // Fuel type 9
    ptr++; strncpy(ptr->fueltype, "C6", 3); ptr->a = 30.0; ptr->b = 0.0800; ptr->c = 3.0;
    ptr->q = 0.80; ptr->bui0 = 62; ptr->cbh = 7; ptr->cfl = 1.80;
    // Fuel type 10
    ptr++; strncpy(ptr->fueltype, "C7", 3); ptr->a = 45.0; ptr->b = 0.0305; ptr->c = 2.0;
    ptr->q = 0.85; ptr->bui0 = 106; ptr->cbh = 10; ptr->cfl = 0.50;
    // Fuel type 11
    ptr++; strncpy(ptr->fueltype, "D1", 3); ptr->a = 30.0; ptr->b = 0.0232; ptr->c = 1.6;
    ptr->q = 0.90; ptr->bui0 = 32; ptr->cbh = 0; ptr->cfl = 0.0;
    // Fuel type 12
    ptr++; strncpy(ptr->fueltype, "S1", 3); ptr->a = 75.0; ptr->b = 0.0297; ptr->c = 1.3;
    ptr->q = 0.75; ptr->bui0 = 38; ptr->cbh = 0; ptr->cfl = 0.0;
    // Fuel type 13
    ptr++; strncpy(ptr->fueltype, "S2", 3); ptr->a = 40.0; ptr->b = 0.0438; ptr->c = 1.7;
    ptr->q = 0.75; ptr->bui0 = 63; ptr->cbh = 0; ptr->cfl = 0.0;
    // Fuel type 14
    ptr++; strncpy(ptr->fueltype, "S3", 3); ptr->a = 55.0; ptr->b = 0.0829; ptr->c = 3.2;
    ptr->q = 0.75; ptr->bui0 = 31; ptr->cbh = 0; ptr->cfl = 0.0;
    // Fuel type 15
    ptr++; strncpy(ptr->fueltype, "O1a", 3); ptr->a = 190.0; ptr->b = 0.0310; ptr->c = 1.40;
    ptr->q = 1.000; ptr->bui0 = 1; ptr->cbh = 0; ptr->cfl = 0.0;
    // Fuel type 16
    ptr++; strncpy(ptr->fueltype, "O1b", 3); ptr->a = 250.0; ptr->b = 0.0350; ptr->c = 1.7;
    ptr->q = 1.000; ptr->bui0 = 1; ptr->cbh = 0; ptr->cfl = 0.0;
    // Fuel type 17
    ptr++; strncpy(ptr->fueltype, "D2", 3); ptr->a = 6.0; ptr->b = 0.0232; ptr->c = 1.6;
    ptr->q = 0.90; ptr->bui0 = 32; ptr->cbh = 0; ptr->cfl = 0.0;
}

// Function to get fuel coefficients for a specific fuel type
fuel_coefs* get_fuel_coefs(const std::string &fueltype, fuel_coefs *fuel_array, int num_fuel_types) {
    for (int i = 0; i < num_fuel_types; ++i) {
        if (fueltype == std::string(fuel_array[i].fueltype)) {
            return &fuel_array[i];
        }
    }
    return nullptr; // Return nullptr if fuel type is not found
}

// Function to parse input CSV
bool parse_input_csv(const std::string &filename, inputs &inp, fuel_coefs *fuel_array, int num_fuel_types) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Unable to open input file: " << filename << std::endl;
        return false;
    }

    std::string line;
    std::getline(file, line); // Skip header line
    if (!std::getline(file, line)) {
        std::cerr << "No data found in file." << std::endl;
        return false;
    }

    std::stringstream ss(line);
    std::string value;

    std::getline(ss, value, ','); 
    std::string trimmed_fueltype = trim(value);  // Make sure trim is a C++98 function
    std::strncpy(inp.fueltype, trimmed_fueltype.c_str(), sizeof(inp.fueltype) - 1);
    inp.fueltype[sizeof(inp.fueltype) - 1] = '\0'; // Ensure null-termination

    std::getline(ss, value, ','); std::getline(ss, value, ','); // Skip columns
    std::getline(ss, value, ','); // Latitude
    std::getline(ss, value, ','); // Longitude
    std::getline(ss, value, ','); // Elevation
    inp.ffmc = safe_stof(value);
    std::getline(ss, value, ','); inp.ws = safe_stof(value);
    std::getline(ss, value, ','); inp.waz = safe_stof(value);
    std::getline(ss, value, ','); inp.bui = safe_stof(value);
    std::getline(ss, value, ','); inp.ps = safe_stof(value);
    std::getline(ss, value, ','); inp.saz = safe_stof(value);
    std::getline(ss, value, ','); inp.pc = safe_stoi(value);
    std::getline(ss, value, ','); inp.pdf = safe_stoi(value);
    std::getline(ss, value, ','); // Skip gfl
    std::getline(ss, value, ','); inp.cur = safe_stof(value);
    std::getline(ss, value, ','); // Skip time
    std::getline(ss, value, ','); // Skip pattern

    fuel_coefs* ptr = get_fuel_coefs(inp.fueltype, fuel_array, num_fuel_types);
    if (ptr == nullptr) {
        std::cerr << "Fuel type not found: " << inp.fueltype << std::endl;
        return false;
    }

    file.close();
    return true;
}

// Function to write output CSV
void write_output_csv(const std::string &filename, float head_ros, float back_ros, float flank_ros) {
    std::ofstream file(filename);
    if (file.is_open()) {
        file << "Head ROS,Back ROS,Flank ROS\n";
        file << head_ros << "," << back_ros << "," << flank_ros << "\n";
        file.close();
    } else {
        std::cerr << "Unable to open output file for writing." << std::endl;
    }
}


int main() {
    inputs inp;
    fuel_coefs fuel_array[18];
    main_outs at;

    // Set up fuel coefficients
    setup_const(fuel_array);

    // Parse input data from CSV
    std::cout << "DEBUGGING : PARSING INPUT CSV!" << std::endl;
    if (!parse_input_csv("/Users/minho/Documents/GitHub/Cell2FireML/notebooks/DataGenerated_Small.csv", inp, fuel_array, 18)) {
        return 1;
    }
    std::cout << "DEBUGGING : DONE PARSING INPUT CSV!" << std::endl;

    // Compute rates of spread
    float head_ros = rate_of_spread(&inp, get_fuel_coefs(inp.fueltype, fuel_array, 18), &at);
    float back_ros = backfire_ros(&inp, get_fuel_coefs(inp.fueltype, fuel_array, 18), &at, at.isi);
    float flank_ros = flankfire_ros(head_ros, back_ros, inp.ws);

    // Save results to CSV file
    write_output_csv("/Users/minho/Documents/GitHub/Cell2FireML/notebooks/DataGenerated_Small_ros.csv", head_ros, back_ros, flank_ros);

    return 0;
}
