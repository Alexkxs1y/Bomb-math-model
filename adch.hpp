#ifndef ADCH_MY
#define ADCH_MY
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <string>

class ADCH{
public:

	ADCH();

	void init(std::string _cx, std::string _cy_alpha, std::string _mx_vr, std::string _mx_wx, std::string _mz_alpha, std::string _mz_wz);

	std::vector<double> alpha_intervals;
	std::vector<double> M_intervals;

	std::vector<double> cx_tabl; // Таблица ад коэф
	std::vector<double> cy_alpha_tabl; // Таблица ад коэф
	std::vector<double> mx_vr_tabl; // Таблица ад коэф
	std::vector<double> mx_wx_tabl; // Таблица ад коэф
	std::vector<double> mz_alpha_tabl; // Таблица ад коэф
	std::vector<double> mz_wz_tabl; // Таблица ад коэф

	double get_cx(double &M, double &alpha);
	double get_cy_alpha(double &M, double &alpha);
	double get_mx_vr(double &M, double &alpha);
	double get_mx_wx(double &M, double &alpha);
	double get_mz_alpha(double &M, double &alpha);
	double get_mz_wz(double &M, double &alpha);
};

#endif