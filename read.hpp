#ifndef READ_MY
#define READ_MY

#include <iostream>
#include <fstream>
#include <string>

void coefsFromFile(std::string filename, std::vector<double> &a);
double findCoef(std::vector<double> &coefs, double &M, double &alpha, std::vector<double> &M_intervals, std::vector<double> &alpha_intervals);

#endif