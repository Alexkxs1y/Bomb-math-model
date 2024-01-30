#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;


void coefsFromFile(string filename,  vector<double> &a){
	ifstream fin(filename);
	for (int i = 0; i < 66; i++){
	        fin >> a[i];
	}
	fin.close(); //закрываем файл
}

double findCoef(vector<double> &coefs, double &M, double &alpha, vector<double> &M_intervals, vector<double> &alpha_intervals){
	int M_nonce = -1;
	int alpha_nonce = -1;
	
	for(int i = 0; i < 11; i++){
		if(M < M_intervals[i]){
			M_nonce = i;
			break;
		}
	}

	for(int i = 0; i < 6; i++){
		if(abs(alpha) < alpha_intervals[i]){
			alpha_nonce = i;
			break;
		}
	}


	if((M_nonce == -1) && (alpha_nonce == -1)){
		return coefs[65];
	}

	if(M_nonce == -1){
		return coefs[11 * (alpha_nonce) - 1] + (coefs[11 * (alpha_nonce + 1) - 1] - coefs[11 * (alpha_nonce) - 1] ) / (alpha_intervals[alpha_nonce] - alpha_intervals[alpha_nonce - 1]) * (alpha - alpha_intervals[alpha_nonce - 1]);
	}

	if(alpha_nonce == -1){
		return coefs[54 + M_nonce] + ( coefs[55 + M_nonce] - coefs[54 + M_nonce] ) / ( M_intervals[M_nonce] - M_intervals[M_nonce - 1] ) * (M - M_intervals[M_nonce - 1]);
	}

	double coef_M0 = coefs[11 * (alpha_nonce - 1) + M_nonce - 1] + ( coefs[11 * (alpha_nonce) + M_nonce - 1] - coefs[11 * (alpha_nonce - 1) + M_nonce - 1] ) /
						( alpha_intervals[alpha_nonce] - alpha_intervals[alpha_nonce - 1] ) *
						( abs(alpha) - alpha_intervals[alpha_nonce - 1]);

	double coef_M1 = coefs[11 * (alpha_nonce - 1) + M_nonce] + ( coefs[11 * (alpha_nonce) + M_nonce] - coefs[11 * (alpha_nonce - 1) + M_nonce] ) /
						( alpha_intervals[alpha_nonce] - alpha_intervals[alpha_nonce - 1] ) *
						( abs(alpha) - alpha_intervals[alpha_nonce - 1]);

	return coef_M0 + (coef_M1 - coef_M0) / (M_intervals[M_nonce] - M_intervals[M_nonce - 1]) * (M - M_intervals[M_nonce - 1]);
}