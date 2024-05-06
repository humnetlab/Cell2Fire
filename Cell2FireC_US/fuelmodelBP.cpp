// Replace the content of this file for each fire spread model (e.g., FBP, Behave, Kitral)
// Helper importation
#include "fuelmodelBP.h"

// Model
#include "rfHROS.h"
#include "rfBROS.h"
#include "rfFROS.h"

// Basic importations
#include <iostream>
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <unordered_map>
#include <string>
#include <vector>
#include <chrono>

// Global Variables for F101
std::unordered_map<int, std::vector<float>> variables;

// Initialize values for FType 101 (sample)
void initialize_var(int nftype)
{
    if (nftype == 101) {
        int F101 = 101;
        std::vector<float> v_101;
        v_101.push_back(0.22);
        v_101.push_back(0);
        v_101.push_back(0);
        v_101.push_back(0.67);
        v_101.push_back(0);
        v_101.push_back(7218);
        v_101.push_back(6562);
        v_101.push_back(0);
        v_101.push_back(12);
        v_101.push_back(15);
        variables.insert(std::make_pair(F101, v_101));
    }
    if (nftype == 102) {
        int F102 = 102;
        std::vector<float> v_102;
        v_102.push_back(0.22);
        v_102.push_back(0);
        v_102.push_back(0);
        v_102.push_back(2.25);
        v_102.push_back(0);
        v_102.push_back(6562);
        v_102.push_back(5906);
        v_102.push_back(0);
        v_102.push_back(30);
        v_102.push_back(15);
        variables.insert(std::make_pair(F102, v_102));
    }
    if (nftype == 103) {
        int F103 = 103;
        std::vector<float> v_103;
        v_103.push_back(0.22);
        v_103.push_back(0.9);
        v_103.push_back(0);
        v_103.push_back(3.37);
        v_103.push_back(0);
        v_103.push_back(4921);
        v_103.push_back(4265);
        v_103.push_back(0);
        v_103.push_back(61);
        v_103.push_back(30);
        variables.insert(std::make_pair(F103, v_103));
    }	
    if (nftype == 104) {
        int F104 = 104;
        std::vector<float> v_104;
        v_104.push_back(0.56);
        v_104.push_back(0);
        v_104.push_back(0);
        v_104.push_back(4.27);
        v_104.push_back(0);
        v_104.push_back(6562);
        v_104.push_back(5906);
        v_104.push_back(0);
        v_104.push_back(61);
        v_104.push_back(15);
        variables.insert(std::make_pair(F104, v_104));
    }
    if (nftype == 105) {
        int F105 = 105;
        std::vector<float> v_105;
        v_105.push_back(0.9);
        v_105.push_back(0);
        v_105.push_back(0);
        v_105.push_back(5.62);
        v_105.push_back(0);
        v_105.push_back(5906);
        v_105.push_back(5249);
        v_105.push_back(0);
        v_105.push_back(46);
        v_105.push_back(40);
        variables.insert(std::make_pair(F105, v_105));
    }
    if (nftype == 106) {
        int F106 = 106;
        std::vector<float> v_106;
        v_106.push_back(0.22);
        v_106.push_back(0);
        v_106.push_back(0);
        v_106.push_back(7.64);
        v_106.push_back(0);
        v_106.push_back(7218);
        v_106.push_back(6562);
        v_106.push_back(0);
        v_106.push_back(46);
        v_106.push_back(40);
        variables.insert(std::make_pair(F106, v_106));
    }	
    if (nftype == 107) {
        int F107 = 107;
        std::vector<float> v_107;
        v_107.push_back(2.25);
        v_107.push_back(0);
        v_107.push_back(0);
        v_107.push_back(12.13);
        v_107.push_back(0);
        v_107.push_back(6562);
        v_107.push_back(5906);
        v_107.push_back(0);
        v_107.push_back(91);
        v_107.push_back(15);
        variables.insert(std::make_pair(F107, v_107));
    }	
        if (nftype == 108) {
        int F108 = 108;
        std::vector<float> v_108;
        v_108.push_back(1.12);
        v_108.push_back(2.25);
        v_108.push_back(0);
        v_108.push_back(16.4);
        v_108.push_back(0);
        v_108.push_back(4921);
        v_108.push_back(4265);
        v_108.push_back(0);
        v_108.push_back(122);
        v_108.push_back(30);
        variables.insert(std::make_pair(F108, v_108));
    }	
    if (nftype == 109) {
        int F109 = 109;
        std::vector<float> v_109;
        v_109.push_back(2.25);
        v_109.push_back(2.25);
        v_109.push_back(0);
        v_109.push_back(20.22);
        v_109.push_back(0);
        v_109.push_back(5906);
        v_109.push_back(5249);
        v_109.push_back(0);
        v_109.push_back(152);
        v_109.push_back(40);
        variables.insert(std::make_pair(F109, v_109));
    }	
    if (nftype == 121) {
        int F121 = 121;
        std::vector<float> v_121;
        v_121.push_back(0.45);
        v_121.push_back(0);
        v_121.push_back(0);
        v_121.push_back(1.12);
        v_121.push_back(1.46);
        v_121.push_back(6562);
        v_121.push_back(5906);
        v_121.push_back(5906);
        v_121.push_back(27);
        v_121.push_back(15);
        variables.insert(std::make_pair(F121, v_121));
    }	
    if (nftype == 122) {
        int F122 = 122;
        std::vector<float> v_122;
        v_122.push_back(1.12);
        v_122.push_back(1.12);
        v_122.push_back(0);
        v_122.push_back(1.35);
        v_122.push_back(2.25);
        v_122.push_back(6562);
        v_122.push_back(5906);
        v_122.push_back(5906);
        v_122.push_back(46);
        v_122.push_back(15);
        variables.insert(std::make_pair(F122, v_122));
    }	
    if (nftype == 123) {
        int F123 = 123;
        std::vector<float> v_123;
        v_123.push_back(0.67);
        v_123.push_back(0.56);
        v_123.push_back(0);
        v_123.push_back(3.26);
        v_123.push_back(2.81);
        v_123.push_back(5906);
        v_123.push_back(5249);
        v_123.push_back(5249);
        v_123.push_back(55);
        v_123.push_back(40);
        variables.insert(std::make_pair(F123, v_123));
    }
    if (nftype == 124) {
        int F124 = 124;
        std::vector<float> v_124;
        v_124.push_back(4.27);
        v_124.push_back(0.67);
        v_124.push_back(0.22);
        v_124.push_back(7.64);
        v_124.push_back(15.96);
        v_124.push_back(5906);
        v_124.push_back(5249);
        v_124.push_back(5249);
        v_124.push_back(64);
        v_124.push_back(40);
        variables.insert(std::make_pair(F124, v_124));
    }	
    if (nftype == 141) {
        int F141 = 141;
        std::vector<float> v_141;
        v_141.push_back(0.56);
        v_141.push_back(0.56);
        v_141.push_back(0);
        v_141.push_back(0.34);
        v_141.push_back(2.92);
        v_141.push_back(6562);
        v_141.push_back(5906);
        v_141.push_back(5249);
        v_141.push_back(30);
        v_141.push_back(15);
        variables.insert(std::make_pair(F141, v_141));
    }
    if (nftype == 142) {
        int F142 = 142;
        std::vector<float> v_142;
        v_142.push_back(3.03);
        v_142.push_back(5.39);
        v_142.push_back(1.69);
        v_142.push_back(0);
        v_142.push_back(8.65);
        v_142.push_back(6562);
        v_142.push_back(0);
        v_142.push_back(5249);
        v_142.push_back(30);
        v_142.push_back(15);
        variables.insert(std::make_pair(F142, v_142));
    }
    if (nftype == 143) {
        int F143 = 143;
        std::vector<float> v_143;
        v_143.push_back(1.01);
        v_143.push_back(6.74);
        v_143.push_back(0);
        v_143.push_back(0);
        v_143.push_back(13.93);
        v_143.push_back(5249);
        v_143.push_back(0);
        v_143.push_back(4593);
        v_143.push_back(73);
        v_143.push_back(40);
        variables.insert(std::make_pair(F143, v_143));
    }
    if (nftype == 144) {
        int F144 = 144;
        std::vector<float> v_144;
        v_144.push_back(1.91);
        v_144.push_back(2.58);
        v_144.push_back(0.45);
        v_144.push_back(0);
        v_144.push_back(5.73);
        v_144.push_back(6562);
        v_144.push_back(5906);
        v_144.push_back(5249);
        v_144.push_back(91);
        v_144.push_back(30);
        variables.insert(std::make_pair(F144, v_144));
    }     
    if (nftype == 145) {
        int F145 = 145;
        std::vector<float> v_145;
        v_145.push_back(8.09);
        v_145.push_back(4.72);
        v_145.push_back(0);
        v_145.push_back(0);
        v_145.push_back(6.52);
        v_145.push_back(2461);
        v_145.push_back(0);
        v_145.push_back(5249);
        v_145.push_back(183);
        v_145.push_back(15);
        variables.insert(std::make_pair(F145, v_145));
    }     
    if (nftype == 146) {
        int F146 = 146;
        std::vector<float> v_146;
        v_146.push_back(6.52);
        v_146.push_back(3.26);
        v_146.push_back(0);
        v_146.push_back(0);
        v_146.push_back(3.15);
        v_146.push_back(2461);
        v_146.push_back(0);
        v_146.push_back(5249);
        v_146.push_back(61);
        v_146.push_back(30);
        variables.insert(std::make_pair(F146, v_146));
    }     
    if (nftype == 147) {
        int F147 = 147;
        std::vector<float> v_147;
        v_147.push_back(7.87);
        v_147.push_back(11.91);
        v_147.push_back(4.94);
        v_147.push_back(0);
        v_147.push_back(7.64);
        v_147.push_back(2461);
        v_147.push_back(0);
        v_147.push_back(5249);
        v_147.push_back(183);
        v_147.push_back(15);
        variables.insert(std::make_pair(F147, v_147));
    }     
    if (nftype == 148) {
        int F148 = 148;
        std::vector<float> v_148;
        v_148.push_back(4.61);
        v_148.push_back(7.64);
        v_148.push_back(1.91);
        v_148.push_back(0);
        v_148.push_back(9.78);
        v_148.push_back(2461);
        v_148.push_back(0);
        v_148.push_back(5249);
        v_148.push_back(91);
        v_148.push_back(40);
        variables.insert(std::make_pair(F148, v_148));
    }     
    if (nftype == 149) {
        int F149 = 149;
        std::vector<float> v_149;
        v_149.push_back(10.11);
        v_149.push_back(5.51);
        v_149.push_back(0);
        v_149.push_back(3.48);
        v_149.push_back(15.73);
        v_149.push_back(2461);
        v_149.push_back(5906);
        v_149.push_back(4921);
        v_149.push_back(134);
        v_149.push_back(40);
        variables.insert(std::make_pair(F149, v_149));
    }
    if (nftype == 144) {
        int F144 = 144;
        std::vector<float> v_144;
        v_144.push_back(1.91);
        v_144.push_back(2.58);
        v_144.push_back(0.45);
        v_144.push_back(0);
        v_144.push_back(5.73);
        v_144.push_back(6562);
        v_144.push_back(5906);
        v_144.push_back(5249);
        v_144.push_back(91);
        v_144.push_back(30);
        variables.insert(std::make_pair(F144, v_144));
    }           
    if (nftype == 161) {
        int F161 = 161;
        std::vector<float> v_161;
        v_161.push_back(0.45);
        v_161.push_back(2.02);
        v_161.push_back(3.37);
        v_161.push_back(0.45);
        v_161.push_back(2.02);
        v_161.push_back(6562);
        v_161.push_back(5906);
        v_161.push_back(5249);
        v_161.push_back(18);
        v_161.push_back(20);
        variables.insert(std::make_pair(F161, v_161));
    }
    if (nftype == 162) {
        int F162 = 162;
        std::vector<float> v_162;
        v_162.push_back(2.13);
        v_162.push_back(4.04);
        v_162.push_back(2.81);
        v_162.push_back(0);
        v_162.push_back(0.45);
        v_162.push_back(6562);
        v_162.push_back(0);
        v_162.push_back(5249);
        v_162.push_back(30);
        v_162.push_back(30);
        variables.insert(std::make_pair(F162, v_162));
    }
    if (nftype == 163) {
        int F163 = 163;
        std::vector<float> v_163;
        v_163.push_back(2.47);
        v_163.push_back(0.34);
        v_163.push_back(0.56);
        v_163.push_back(1.46);
        v_163.push_back(2.47);
        v_163.push_back(5906);
        v_163.push_back(5249);
        v_163.push_back(4593);
        v_163.push_back(40);
        v_163.push_back(30);
        variables.insert(std::make_pair(F163, v_163));
    }
    if (nftype == 164) {
        int F164 = 164;
        std::vector<float> v_164;
        v_164.push_back(10.11);
        v_164.push_back(0);
        v_164.push_back(0);
        v_164.push_back(0);
        v_164.push_back(4.49);
        v_164.push_back(7546);
        v_164.push_back(0);
        v_164.push_back(6562);
        v_164.push_back(15);
        v_164.push_back(12);
        variables.insert(std::make_pair(F164, v_164));
    }
    if (nftype == 165) {
        int F165 = 165;
        std::vector<float> v_165;
        v_165.push_back(8.99);
        v_165.push_back(8.99);
        v_165.push_back(6.74);
        v_165.push_back(0);
        v_165.push_back(6.74);
        v_165.push_back(4921);
        v_165.push_back(0);
        v_165.push_back(2461);
        v_165.push_back(30);
        v_165.push_back(25);
        variables.insert(std::make_pair(F165, v_165));
    }
    if (nftype == 181) {
        int F181 = 181;
        std::vector<float> v_181;
        v_181.push_back(2.25);
        v_181.push_back(4.94);
        v_181.push_back(8.09);
        v_181.push_back(0);
        v_181.push_back(0);
        v_181.push_back(6562);
        v_181.push_back(0);
        v_181.push_back(0);
        v_181.push_back(6);
        v_181.push_back(30);
        variables.insert(std::make_pair(F181, v_181));
    }
    if (nftype == 182) {
        int F182 = 182;
        std::vector<float> v_182;
        v_182.push_back(3.15);
        v_182.push_back(5.17);
        v_182.push_back(4.94);
        v_182.push_back(0);
        v_182.push_back(0);
        v_182.push_back(6562);
        v_182.push_back(0);
        v_182.push_back(0);
        v_182.push_back(6);
        v_182.push_back(25);
        variables.insert(std::make_pair(F182, v_182));
    }
    if (nftype == 183) {
        int F183 = 183;
        std::vector<float> v_183;
        v_183.push_back(1.12);
        v_183.push_back(4.94);
        v_183.push_back(6.29);
        v_183.push_back(0);
        v_183.push_back(0);
        v_183.push_back(6562);
        v_183.push_back(0);
        v_183.push_back(0);
        v_183.push_back(9);
        v_183.push_back(20);
        variables.insert(std::make_pair(F183, v_183));
    }
    if (nftype == 184) {
        int F184 = 184;
        std::vector<float> v_184;
        v_184.push_back(1.12);
        v_184.push_back(3.37);
        v_184.push_back(9.44);
        v_184.push_back(0);
        v_184.push_back(0);
        v_184.push_back(6562);
        v_184.push_back(0);
        v_184.push_back(0);
        v_184.push_back(12);
        v_184.push_back(25);
        variables.insert(std::make_pair(F184, v_184));
    }
    if (nftype == 185) {
        int F185 = 185;
        std::vector<float> v_185;
        v_185.push_back(2.58);
        v_185.push_back(5.62);
        v_185.push_back(9.89);
        v_185.push_back(0);
        v_185.push_back(0);
        v_185.push_back(6562);
        v_185.push_back(0);
        v_185.push_back(5249);
        v_185.push_back(18);
        v_185.push_back(25);
        variables.insert(std::make_pair(F185, v_185));
    }
    if (nftype == 186) {
        int F186 = 186;
        std::vector<float> v_186;
        v_186.push_back(5.39);
        v_186.push_back(2.7);
        v_186.push_back(2.7);
        v_186.push_back(0);
        v_186.push_back(0);
        v_186.push_back(6562);
        v_186.push_back(0);
        v_186.push_back(0);
        v_186.push_back(9);
        v_186.push_back(25);
        variables.insert(std::make_pair(F186, v_186));
    }
    if (nftype == 187) {
        int F187 = 187;
        std::vector<float> v_187;
        v_187.push_back(0.67);
        v_187.push_back(3.15);
        v_187.push_back(18.2);
        v_187.push_back(0);
        v_187.push_back(0);
        v_187.push_back(6562);
        v_187.push_back(0);
        v_187.push_back(0);
        v_187.push_back(12);
        v_187.push_back(25);
        variables.insert(std::make_pair(F187, v_187));
    }
    if (nftype == 188) {
        int F188 = 188;
        std::vector<float> v_188;
        v_188.push_back(13.03);
        v_188.push_back(3.15);
        v_188.push_back(2.47);
        v_188.push_back(0);
        v_188.push_back(0);
        v_188.push_back(5906);
        v_188.push_back(0);
        v_188.push_back(0);
        v_188.push_back(9);
        v_188.push_back(35);
        variables.insert(std::make_pair(F188, v_188));
    }
    if (nftype == 189) {
        int F189 = 189;
        std::vector<float> v_189;
        v_189.push_back(14.94);
        v_189.push_back(7.42);
        v_189.push_back(9.33);
        v_189.push_back(0);
        v_189.push_back(0);
        v_189.push_back(5906);
        v_189.push_back(0);
        v_189.push_back(5249);
        v_189.push_back(18);
        v_189.push_back(35);
        variables.insert(std::make_pair(F189, v_189));
    }
    if (nftype == 201) {
        int F201 = 201;
        std::vector<float> v_201;
        v_201.push_back(3.37);
        v_201.push_back(6.74);
        v_201.push_back(24.72);
        v_201.push_back(0);
        v_201.push_back(0);
        v_201.push_back(6562);
        v_201.push_back(0);
        v_201.push_back(0);
        v_201.push_back(30);
        v_201.push_back(25);
        variables.insert(std::make_pair(F201, v_201));
    }
    if (nftype == 202) {
        int F202 = 202;
        std::vector<float> v_202;
        v_202.push_back(10.11);
        v_202.push_back(9.55);
        v_202.push_back(8.99);
        v_202.push_back(0);
        v_202.push_back(0);
        v_202.push_back(6562);
        v_202.push_back(0);
        v_202.push_back(0);
        v_202.push_back(30);
        v_202.push_back(25);
        variables.insert(std::make_pair(F202, v_202));
    }
    if (nftype == 203) {
        int F203 = 203;
        std::vector<float> v_203;
        v_203.push_back(12.36);
        v_203.push_back(6.18);
        v_203.push_back(6.74);
        v_203.push_back(0);
        v_203.push_back(0);
        v_203.push_back(6562);
        v_203.push_back(0);
        v_203.push_back(0);
        v_203.push_back(37);
        v_203.push_back(25);
        variables.insert(std::make_pair(F203, v_203));
    }
    if (nftype == 204) {
        int F204 = 204;
        std::vector<float> v_204;
        v_204.push_back(11.8);
        v_204.push_back(7.87);
        v_204.push_back(11.8);
        v_204.push_back(0);
        v_204.push_back(0);
        v_204.push_back(6562);
        v_204.push_back(0);
        v_204.push_back(0);
        v_204.push_back(82);
        v_204.push_back(25);
        variables.insert(std::make_pair(F204, v_204));
    }
    if (nftype == 91) {
        int F91 = 91;
        std::vector<float> v_91;
        v_91.push_back(0);
        v_91.push_back(0);
        v_91.push_back(0);
        v_91.push_back(0);
        v_91.push_back(0);
        v_91.push_back(0);
        v_91.push_back(0);
        v_91.push_back(0);
        v_91.push_back(0);
        v_91.push_back(0);
        variables.insert(std::make_pair(F91, v_91));
    }
    if (nftype == 92) {
        int F92 = 92;
        std::vector<float> v_92;
        v_92.push_back(0);
        v_92.push_back(0);
        v_92.push_back(0);
        v_92.push_back(0);
        v_92.push_back(0);
        v_92.push_back(0);
        v_92.push_back(0);
        v_92.push_back(0);
        v_92.push_back(0);
        v_92.push_back(0);
        variables.insert(std::make_pair(F92, v_92));
    }
    if (nftype == 93) {
        int F93 = 93;
        std::vector<float> v_93;
        v_93.push_back(0);
        v_93.push_back(0);
        v_93.push_back(0);
        v_93.push_back(0);
        v_93.push_back(0);
        v_93.push_back(0);
        v_93.push_back(0);
        v_93.push_back(0);
        v_93.push_back(0);
        v_93.push_back(0);
        variables.insert(std::make_pair(F93, v_93));
    }
    if (nftype == 98) {
        int F98 = 98;
        std::vector<float> v_98;
        v_98.push_back(0);
        v_98.push_back(0);
        v_98.push_back(0);
        v_98.push_back(0);
        v_98.push_back(0);
        v_98.push_back(0);
        v_98.push_back(0);
        v_98.push_back(0);
        v_98.push_back(0);
        v_98.push_back(0);
        variables.insert(std::make_pair(F98, v_98));
    }
    if (nftype == 99) {
        int F99 = 99;
        std::vector<float> v_99;
        v_99.push_back(0);
        v_99.push_back(0);
        v_99.push_back(0);
        v_99.push_back(0);
        v_99.push_back(0);
        v_99.push_back(0);
        v_99.push_back(0);
        v_99.push_back(0);
        v_99.push_back(0);
        v_99.push_back(0);
        variables.insert(std::make_pair(F99, v_99));
    }
    if (nftype == 0) {
        int F0 = 0;
        std::vector<float> v_0;
        v_0.push_back(0);
        v_0.push_back(0);
        v_0.push_back(0);
        v_0.push_back(0);
        v_0.push_back(0);
        v_0.push_back(0);
        v_0.push_back(0);
        v_0.push_back(0);
        v_0.push_back(0);
        v_0.push_back(0);
        variables.insert(std::make_pair(F0, v_0));
    }
}


/* ----------------- Length-to-Breadth --------------------------*/
// Returns length-to-breadth ratio given ws in km/hr (factor transforms to m/s inside
float l_to_b(float ws)
{
    float alpha, beta, factor;
    alpha = 0.2566;
    beta = -0.1548;
    factor = 1000.0 / 3600.0;
    return pow((0.936 * exp(alpha * factor * ws) + 0.461 * exp(beta * factor * ws) - 0.397), 0.45);
}


/* Functions */
// Calculates ROS in three main directions, LB, and elliptical components (a,b,c) fiven input data
// Calls the ML model of preference (RF in this version)
void calculate_BP(inputs *data, 
                  fuel_coefs *ptr,
                  main_outs *at,
                  snd_outs *sec,
                  fire_struc *hptr,
                  fire_struc *fptr,
                  fire_struc *bptr)
{
    /* (1) Input parameters */
    initialize_var(data->nftype);

    // Input variables (8)
    float mc1, mc10, mc100, mcWoody, mcHerb; // MC is traditionally int, but we use float here for precision
    float ws, waz, slope;
    int nftype;

    // Grab data from inputs
    nftype = data->nftype;
    ws = data->ws;
    waz = data->waz;
    slope = data->slope;
    mc1 = data->mc1;
    mc10 = data->mc10;
    mc100 = data->mc100;
    mcHerb = data->mcHerb;
    mcWoody = data->mcWoody;

    // Print outputs (Sanity check)
    
    // std::cout << "---------------------------------------" << std::endl;
    // std::cout << "(BP)Sanity Check for variables" << std::endl;
    // std::cout << "nftype : " << data->nftype << std::endl;
    // std::cout << "ws : " << ws << ", " << "waz : " << waz << std::endl;
    // std::cout << "slope : " << slope << std::endl;
    // std::cout << "mc1 : " << mc1 << ", " << "mc10 : " << mc10 << ", " << "mc100 : " << mc100 << std::endl;
    // std::cout << "mcWoody : " << mcWoody << ", " << "mcHerb : " << mcHerb << std::endl;
    // std::cout << "---------------------------------------" << std::endl;
    

    // Set pointers
    ptr->nftype = data->nftype;
    ptr->ws=ws; 
    ptr->waz=waz;
    ptr->slope=slope; 
    ptr->mc1=mc1;
    ptr->mc10=mc10;
    ptr->mc100=mc100;
    ptr->mcWoody=mcWoody;
    ptr->mcHerb=mcHerb;

    // Create input tensor from sample row
    double array[] = {static_cast<double>(data->nftype),
                      ws,
                      slope,
                      mc1,
                      mc10,
                      mc100,
                      mcWoody,
                      mcHerb, 
                      variables[data->nftype][0],
                      variables[data->nftype][1],
                      variables[data->nftype][2],
                      variables[data->nftype][3],
                      variables[data->nftype][4],
                      variables[data->nftype][5],
                      variables[data->nftype][6],
                      variables[data->nftype][7],
                      variables[data->nftype][8],
                      variables[data->nftype][9]};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /*****************
       Start time
    *****************/    
    auto start = std::chrono::high_resolution_clock::now();

    // Model predict
    double hros_pred = hros(array);
    double bros_pred = bros(array);
    double fros_pred = fros(array);	
    
    /*****************
        End time
    *****************/
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "\n(BP)Predicted : " << hros_pred << " | " << bros_pred << " | " << fros_pred << std::endl;
    
    // Calculate the duration by subtracting start time from end time
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    // Print the duration in nanoseconds
    std::cout << "Time taken: " << duration.count() << " nanoseconds" << std::endl;


    // Create HROS, BROS, FROS to pass
    auto hrss = hros_pred;
    auto brss = bros_pred;
    auto frss = fros_pred;

    // Outputs to pointers (rss and ros, no difference in this version)
    at->rss = hrss;
    hptr->rss = hrss;
    bptr->rss = brss;
    fptr->rss = frss;

    hptr->ros = hptr->rss;
    bptr->ros = bptr->rss; 
    fptr->ros = fptr->rss;

    // Ellipse components from basic function (LB, a, b, c)
    sec->lb = l_to_b(data->ws);
    at->a = (hptr->rss + bptr->rss) / 2.;
    at->b = (hptr->rss + bptr->rss) / (2. * sec->lb); 
    at->c = (hptr->rss - bptr->rss) / 2.; 
    //std::cout << "(BP)lb: " << sec->lb << " | a: " << at->a << " | b: " << at->b << " | c: " << at->c << std::endl;
}
