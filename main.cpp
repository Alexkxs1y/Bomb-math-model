#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include "mina.hpp"

using namespace std;

int main(){
	double m = 20;
	vector<double> J = {0.041, 0.753, 0.753};
	double d = 0.11;
	double l = 0.9;
	double V0 = 441.7;
	double Vx = V0 * cos(50.0 * M_PI / 180.0);
	double Vy = V0 * sin(50.0 * M_PI / 180.0);
	vector<double> stateVectorG = {0, 0.0001, 0, Vx , Vy, 0};
	vector<double> angles = {0, 0, 50.0 * M_PI / 180.0};
	vector<double> speedAngles = {0, 0, 0};
	vector<double> angleSpeed = {0, 0, 0};
	vector<double> targetPosition = {7921.86, 137.61};
	double li = 0.004;
	double tEng = 0.03;
	double engThrust = 4903.33;
	double t = 0;
	double t_last = 0;
	bool corr = false;
	double dt = 0.0001;
	int iter = 0;

	GHOST4401 atm;
	atm.init();

	ADCH adch_start;
	adch_start.init("cx.dat", "cy_alpha.dat", "mx_vr.dat", "mx_wx.dat", "mz_alpha.dat", "mz_wz.dat");

	ADCH adch_corr;
	adch_corr.init("cx_corr.dat", "cy_alpha_corr.dat", "mx_vr_corr.dat", "mx_wx_corr.dat", "mz_alpha_corr.dat", "mz_wz_corr.dat");

	otr mina;
	mina.init(m, J, d, l, li, tEng, engThrust, t, t_last, stateVectorG, angles, speedAngles, angleSpeed, targetPosition, adch_start, atm, corr);

	vector<double> res(3);
	res = mina.get_res(dt, iter);
	iter = 0;
	double t_of_flight = mina.t;

	mina.init(m, J, d, l, li, tEng, engThrust, t, t_last, stateVectorG, angles, speedAngles, angleSpeed, targetPosition, adch_start, atm, corr);

	//Полет до сброса обтекателя
	while(mina.t < t_of_flight - 3.5){
		mina.update(dt, iter);
	}

	//Сброс обтекателя
	double m_corr = 19;
	vector<double> J_corr = {0.0407, 0.576, 0.576};
	for(int i = 0; i < 3; i++){
		angleSpeed[i] = mina.orientationVector[i + 4];
	}
	mina.init(m_corr, J_corr, d, l, li, tEng, engThrust, mina.t, t_last, mina.stateVectorG, mina.angleVector, mina.speedAngles,
			 angleSpeed, mina.Rg, adch_corr, atm, corr);

	//Полёт без обтекателя
	while(mina.t < t_of_flight - 3.0){
		mina.update(dt, iter);
	}

	//Включение ИК
	corr = true;
	for(int i = 0; i < 3; i++){
			angleSpeed[i] = mina.orientationVector[i + 4];
	}
	double maxRange_x = 260;
	double maxRange_z = 260;
	double step_x = 4;
	double step_z = 4;
	double t_corr = mina.t;
	vector<double> stateVectorG_corr = mina.stateVectorG;
	vector<double> angleVector_corr = mina.angleVector;
	vector<double> speedAngles_corr = mina.speedAngles;

	for(int i = 0; i < 131; i++){
		targetPosition[0] = res[0] - maxRange_x + step_x * double(i);
		for(int j = 0; j < 131; j++){
			targetPosition[1] = res[2] - maxRange_z + step_z * double(j);	
			mina.init(m_corr, J_corr, d, l, li, tEng, engThrust, t_corr, t_last, stateVectorG_corr, angleVector_corr, speedAngles_corr,
					 angleSpeed, targetPosition, adch_corr, atm, corr);
			vector<double> res_new(3);
			res_new = mina.get_res(dt);
			if( sqrt( (res_new[0] - targetPosition[0])*(res_new[0] - targetPosition[0]) + (res_new[2] - targetPosition[1])*(res_new[2] - targetPosition[1]) ) < 50 ){
				cout << targetPosition[0] << ' ' << targetPosition[1] << ' ' << res_new[0] << ' ' << res_new[2] << '\n';
			}
		}
	}

	return 1;
}