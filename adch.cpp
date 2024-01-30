#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include "adch.hpp"
#include "read.hpp"

ADCH::ADCH():cx_tabl(66), cy_alpha_tabl(66), mx_vr_tabl(66), mx_wx_tabl(66), mz_alpha_tabl(66), mz_wz_tabl(66){}

void ADCH::init(std::string _cx, std::string _cy_alpha, std::string _mx_vr, std::string _mx_wx, std::string _mz_alpha, std::string _mz_wz){
	coefsFromFile(_cx, cx_tabl);
	coefsFromFile(_cy_alpha, cy_alpha_tabl);
	coefsFromFile(_mx_vr, mx_vr_tabl);
	coefsFromFile(_mx_wx, mx_wx_tabl);
	coefsFromFile(_mz_alpha, mz_alpha_tabl);
	coefsFromFile(_mz_wz, mz_wz_tabl);
	alpha_intervals = {0.0 , 1.0 , 2.0 , 3.0 , 4.0, 5.0};
	M_intervals = {.000, .500, .600, .700, .800, .900, 1.000, 1.100, 1.200, 1.300, 1.400};
}

double ADCH::get_cx(double &M, double &alpha){
	return findCoef(cx_tabl, M, alpha, M_intervals, alpha_intervals);
}

double ADCH::get_cy_alpha(double &M, double &alpha){
	return findCoef(cy_alpha_tabl, M, alpha, M_intervals, alpha_intervals);
}

double ADCH::get_mx_vr(double &M, double &alpha){
	return findCoef(mx_vr_tabl, M, alpha, M_intervals, alpha_intervals);
}

double ADCH::get_mx_wx(double &M, double &alpha){
	return findCoef(mx_wx_tabl, M, alpha, M_intervals, alpha_intervals);
}

double ADCH::get_mz_alpha(double &M, double &alpha){
	return findCoef(mz_alpha_tabl, M, alpha, M_intervals, alpha_intervals);
}

double ADCH::get_mz_wz(double &M, double &alpha){
	return findCoef(mz_wz_tabl, M, alpha, M_intervals, alpha_intervals);
}