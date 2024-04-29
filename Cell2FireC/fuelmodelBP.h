// Helper definition
#ifndef fuelmodelBP
#define fuelmodelBP

// Basic importations
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <unordered_map>
#include <string>
#include <vector>

// Torch (if needed)
#include <torch/torch.h>
#include <torch/script.h> 

/*
    Input, Fire, and output structures
*/

typedef struct
   { char fueltype[4];
     float ws, saz, cur, ps, cbh, ccf, cbd, elev, slope, mc1, mc10, mc100, mcWoody, mcHerb ;
     int waz, nftype ;
     float factor_cbd, factor_ccf, factor_ros10, factor_actv ;
     bool cros, verbose;
   } inputs;
  
  typedef struct
   { float rss, a, b, c, angle; 
     char covertype;
     bool cros;
   } main_outs;

typedef struct
   { char fueltype[4] ;
     float ws, waz, slope, mc1, mc10, mc100, mcWoody, mcHerb;
     int nftype;
   } fuel_coefs;

typedef struct
   { float lb ;
   } snd_outs;

typedef struct
   { float ros, rss, brss, frss, hros, bros, fros;
   } fire_struc;


 /*
    Functions
 */

// Initialize coefficients
void initialize_var(int nftype);
 
// Length-to-Breadth ratio
float l_to_b(float ws);

// Main function to populate spread outputs based on inputs provided from main class
// void calculate_BP(inputs *data, main_outs *at, fire_struc *hptr);
void calculate_BP(inputs *data,
                  fuel_coefs *ptr, 
                  main_outs *at,
                  snd_outs *sec,
                  fire_struc *hptr,
                  fire_struc *fptr,
                  fire_struc *bptr);

#endif