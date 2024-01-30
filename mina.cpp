#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include "mina.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////
//						Конструктор класса						 //
///////////////////////////////////////////////////////////////////

otr::otr():J(3), engNumber(10), stateVectorG(6), angleVector(6), d_stateVectorG(6), orientationVector(7), d_orientationVector(7),
				speedAngles(2), bg(3), bi(3), Sg(3), Rg(3), phi_i(2), g(3), aeroForces(3), aeroForces_G(3), aeroTorque(3), thrustForce(3){}


///////////////////////////////////////////////////////////////////
//						Инициализация класса					 //
///////////////////////////////////////////////////////////////////

void otr::init(	double &_m, vector<double> &_J, double &_d, double &_l, double &_li, double &_tEng, 
				double &_engThrust, double &_t, double &_t_last, vector<double> &_stateVectorG, 
				vector<double> &_angles, vector<double> &_speedAngles, vector<double> &_angleSpeed, vector<double> &_targetPosition,
				ADCH &_adch, GHOST4401 &_atm, bool _corr)
{

	m = _m; J = _J; d = _d; l = _l; li = _li; tEng = _tEng; engThrust = _engThrust; t = _t; t_last = _t_last; stateVectorG = _stateVectorG; angleVector = _angles;
	speedAngles = _speedAngles; atm = _atm; adch = _adch, corr = _corr;

	Rg[0] = _targetPosition[0];
	Rg[1] = 0;
	Rg[2] = _targetPosition[1];
	workEngNumber = -1;

	for(int i = 0; i < 10; i++){
		engNumber[i] = 0;
	}

	for(int i = 0; i < 3; i ++){
		orientationVector[i + 4] = _angleSpeed[i];
	}

	first_init_RodGam();
	init_step_params();
}


///////////////////////////////////////////////////////////////////
//					Родриго из углов Эйлера						 //
///////////////////////////////////////////////////////////////////

void otr::first_init_RodGam(){
	orientationVector[0] = cos(angleVector[1] / 2) * cos(angleVector[0] / 2) * cos(angleVector[2] / 2) - sin(angleVector[1] / 2) * sin(angleVector[2] / 2) * sin(angleVector[0] / 2);
	orientationVector[1] = sin(angleVector[1] / 2) * cos(angleVector[0] / 2) * sin(angleVector[2] / 2) + cos(angleVector[1] / 2) * cos(angleVector[2] / 2) * sin(angleVector[0] / 2);
	orientationVector[2] = sin(angleVector[1] / 2) * cos(angleVector[2] / 2) * cos(angleVector[0] / 2) + cos(angleVector[1] / 2) * sin(angleVector[2] / 2) * sin(angleVector[0] / 2);
	orientationVector[3] = cos(angleVector[1] / 2) * sin(angleVector[2] / 2) * cos(angleVector[0] / 2) - sin(angleVector[1] / 2) * cos(angleVector[2] / 2) * sin(angleVector[0] / 2);
	double norm = sqrt(orientationVector[0]*orientationVector[0] + orientationVector[1]*orientationVector[1] +
						orientationVector[2]*orientationVector[2] + orientationVector[3]*orientationVector[3]);
	for(int i = 0; i < 4; i ++){
		orientationVector[i] /= norm;
	}
}



///////////////////////////////////////////////////////////////////
//					Дополнительные параметры					 //
///////////////////////////////////////////////////////////////////

void otr::init_extra_params(){
	Vref = sqrt( (stateVectorG[3])*(stateVectorG[3]) + (stateVectorG[4])*(stateVectorG[4]) + 
				(stateVectorG[5] )*(stateVectorG[5])  );
	M = Vref / atm.get_a(stateVectorG[1]);
	q = Vref * Vref * atm.get_density(stateVectorG[1]) * 0.5;
	qs = q*M_PI*d*d/4;
	qsl = qs*l;
}

/////////////////////////////////////////////////////////////////////////
//				Определение углов и угловых скоростей				   //
/////////////////////////////////////////////////////////////////////////

void otr::init_angleVector(){
	angleVector[0] = atan((2 * (orientationVector[0] * orientationVector[1] - orientationVector[2] * orientationVector[3])) / 
								(2 * (orientationVector[0] * orientationVector[0] + orientationVector[2] * orientationVector[2]) - 1));
	angleVector[1] = atan((2 * (orientationVector[0] * orientationVector[2] - orientationVector[3] * orientationVector[1])) / 
								(2 * (orientationVector[0] * orientationVector[0] + orientationVector[1] * orientationVector[1]) - 1));
	angleVector[2] = asin(2 * (orientationVector[0] * orientationVector[3] + orientationVector[1] * orientationVector[2]));
	
	angleVector[3] = orientationVector[4] - (tan(angleVector[2]) * (orientationVector[5] * cos(angleVector[0]) - orientationVector[6] * sin(angleVector[0])));
	angleVector[4] = ((1 / cos(angleVector[2])) * (orientationVector[5] * cos(angleVector[0]) - orientationVector[6] * sin(angleVector[0])));
	angleVector[5] = orientationVector[5] * sin(angleVector[0]) + orientationVector[6] * cos(angleVector[0]);
}

void otr::init_speedAngles(){
	double a11 = 2 * (orientationVector[0] * orientationVector[0] + orientationVector[1] * orientationVector[1]) - 1;
	double a12 = 2 * (- orientationVector[0] * orientationVector[3] + orientationVector[1] * orientationVector[2]);
	double a13 = 2 * (orientationVector[0] * orientationVector[2] + orientationVector[1] * orientationVector[3]);

	double a21 = 2 * (orientationVector[0] * orientationVector[3] + orientationVector[1] * orientationVector[2]);
	double a22 = 2 * (orientationVector[0] * orientationVector[0] + orientationVector[2] * orientationVector[2]) - 1;
	double a23 = 2 * (- orientationVector[0] * orientationVector[1] + orientationVector[2] * orientationVector[3]);

	double a31 = 2 * (- orientationVector[0] * orientationVector[2] + orientationVector[1] * orientationVector[3]);
	double a32 = 2 * (orientationVector[0] * orientationVector[1] + orientationVector[2] * orientationVector[3]);
	double a33 = 2 * (orientationVector[0] * orientationVector[0] + orientationVector[3] * orientationVector[3]) - 1;

	double Vx = a11 * stateVectorG[3] + a21 * stateVectorG[4] + a31 * stateVectorG[5];
	double Vy = a12 * stateVectorG[3] + a22 * stateVectorG[4] + a32 * stateVectorG[5];
	double Vz = a13 * stateVectorG[3] + a23 * stateVectorG[4] + a33 * stateVectorG[5];

	speedAngles[0] = - atan(Vy/Vx);
	speedAngles[1] = asin(Vz/Vref);
}

/////////////////////////////////////////////////////////////////////////
//					Поиск углов цели на визоре						   //
/////////////////////////////////////////////////////////////////////////

void otr::init_b(){
	Sg[0] = li * cos(angleVector[1] + speedAngles[1]) * cos(angleVector[2] + speedAngles[0]);
	Sg[1] = li * cos(angleVector[1] + speedAngles[1]) * sin(angleVector[2] + speedAngles[0]);
	Sg[2] = li * sin(angleVector[1] + speedAngles[1]);
	for(int i = 0; i < 3; i++){
		bg[i] = Rg[i] - stateVectorG[i] - Sg[i];
	}

	double a11 = 2 * (orientationVector[0] * orientationVector[0] + orientationVector[1] * orientationVector[1]) - 1;
	double a12 = 2 * (- orientationVector[0] * orientationVector[3] + orientationVector[1] * orientationVector[2]);
	double a13 = 2 * (orientationVector[0] * orientationVector[2] + orientationVector[1] * orientationVector[3]);

	double a21 = 2 * (orientationVector[0] * orientationVector[3] + orientationVector[1] * orientationVector[2]);
	double a22 = 2 * (orientationVector[0] * orientationVector[0] + orientationVector[2] * orientationVector[2]) - 1;
	double a23 = 2 * (- orientationVector[0] * orientationVector[1] + orientationVector[2] * orientationVector[3]);

	double a31 = 2 * (- orientationVector[0] * orientationVector[2] + orientationVector[1] * orientationVector[3]);
	double a32 = 2 * (orientationVector[0] * orientationVector[1] + orientationVector[2] * orientationVector[3]);
	double a33 = 2 * (orientationVector[0] * orientationVector[0] + orientationVector[3] * orientationVector[3]) - 1;

	bi[0] = a11 * bg[0] + a21 * bg[1] + a31 * bg[2];
	bi[1] = a12 * bg[0] + a22 * bg[1] + a32 * bg[2];
	bi[2] = a13 * bg[0] + a23 * bg[1] + a33 * bg[2];
}

void otr::init_phi_i(){
	phi_i[0] = atan(sqrt(bi[1]*bi[1] + bi[2]*bi[2]) / bi[0]);
	phi_i[1] = atan2(bi[2],bi[1]);
}


///////////////////////////////////////////////////////////////////
//						Силы и моменты		 					 //
///////////////////////////////////////////////////////////////////

void otr::init_g(){
	g[0] = g[2] = 0;
	g[1] = atm.get_g(stateVectorG[1]);
}

void otr::init_aeroForces(){
	double alpha_deg = 180.0 * speedAngles[0] / M_PI;
	double beta_deg = 180.0 * speedAngles[1] / M_PI;
	aeroForces[0] = - adch.get_cx(M, alpha_deg) * qs;
	aeroForces[1] = adch.get_cy_alpha(M, alpha_deg) * speedAngles[0] * qs;
	aeroForces[2] = - adch.get_cy_alpha(M, beta_deg) * speedAngles[1] * qs ;
}

void otr::init_aeroForces_G(){
	double a11 = 2 * (orientationVector[0] * orientationVector[0] + orientationVector[1] * orientationVector[1]) - 1;
	double a12 = 2 * (- orientationVector[0] * orientationVector[3] + orientationVector[1] * orientationVector[2]);
	double a13 = 2 * (orientationVector[0] * orientationVector[2] + orientationVector[1] * orientationVector[3]);

	double a21 = 2 * (orientationVector[0] * orientationVector[3] + orientationVector[1] * orientationVector[2]);
	double a22 = 2 * (orientationVector[0] * orientationVector[0] + orientationVector[2] * orientationVector[2]) - 1;
	double a23 = 2 * (- orientationVector[0] * orientationVector[1] + orientationVector[2] * orientationVector[3]);

	double a31 = 2 * (- orientationVector[0] * orientationVector[2] + orientationVector[1] * orientationVector[3]);
	double a32 = 2 * (orientationVector[0] * orientationVector[1] + orientationVector[2] * orientationVector[3]);
	double a33 = 2 * (orientationVector[0] * orientationVector[0] + orientationVector[3] * orientationVector[3]) - 1;

	if(corr){
		aeroForces_G[0] = a11 * (aeroForces[0]) + a12 * (aeroForces[1] + thrustForce[1]) + a13 * (aeroForces[2] + thrustForce[2]);
		aeroForces_G[1] = a21 * (aeroForces[0]) + a22 * (aeroForces[1] + thrustForce[1]) + a23 * (aeroForces[2] + thrustForce[2]);
		aeroForces_G[2] = a31 * (aeroForces[0]) + a32 * (aeroForces[1] + thrustForce[1]) + a33 * (aeroForces[2] + thrustForce[2]);
	} else {
		aeroForces_G[0] = a11 * aeroForces[0] + a12 * aeroForces[1] + a13 * aeroForces[2];
		aeroForces_G[1] = a21 * aeroForces[0] + a22 * aeroForces[1] + a23 * aeroForces[2];
		aeroForces_G[2] = a31 * aeroForces[0] + a32 * aeroForces[1] + a33 * aeroForces[2];
	}
}

void otr::init_aeroTorque(){
	double alpha_deg = 180.0 * speedAngles[0] / M_PI;
	double beta_deg = 180.0 * speedAngles[1] / M_PI;
	aeroTorque[0] = (adch.get_mx_vr(M, alpha_deg) + adch.get_mx_wx(M, alpha_deg) * l / Vref * orientationVector[4]) * qsl;
	aeroTorque[1] = (adch.get_mz_alpha(M, beta_deg) * speedAngles[1] + adch.get_mz_wz(M, beta_deg) * l / Vref * orientationVector[5]) * qsl;
	aeroTorque[2] = (adch.get_mz_alpha(M, alpha_deg) * speedAngles[0] + adch.get_mz_wz(M, alpha_deg) * l / Vref * orientationVector[6]) * qsl;
}

void otr::init_forceAndTorque(){
	init_g();
	init_aeroForces();
	init_thrustForce();
	init_aeroForces_G();
	init_aeroTorque();
}

void otr::init_thrustForce(){
	if(t - t_last > tEng){
		thrustForce[0] = 0;
		thrustForce[1] = 0;
		thrustForce[2] = 0;

		double del = (phi_i[1] - orientationVector[4] / 2 * tEng);
		if(del > 2.0 * M_PI){
			del -= 2.0 * M_PI;
		}
		if(del < 0){
			del += 2.0 * M_PI;
		}

		for(int i = 0; i < 10; i++){
			if((engNumber[(i+5)%10] >= 0) && (abs(del - 2.0 * M_PI / 10.0 * double(i)) < 0.04) && (abs(phi_i[0]) > 5.0 * M_PI / 180.0) && ((abs(phi_i[0]) < 13.0 * M_PI / 180.0)) ){
				workEngNumber = i ;
				t_last = t;
				thrustForce[0] = 0;
				thrustForce[1] = engThrust * cos(2.0 * M_PI / 10.0 * double(workEngNumber));
				thrustForce[2] = engThrust * sin(2.0 * M_PI / 10.0 * double(workEngNumber));
				engNumber[(i+5)%10] = -2;
				break;
			}
		}
	}	
}


///////////////////////////////////////////////////////////////////
//						Производные состояния 					 //
///////////////////////////////////////////////////////////////////

void otr::calc_d_stateVectorG(){
	for(int i = 0; i < 3; i++){
		d_stateVectorG[i] = stateVectorG[3 + i];
	}
	d_stateVectorG[3] = aeroForces_G[0] / m;
	d_stateVectorG[4] = aeroForces_G[1] / m - g[1];
	d_stateVectorG[5] = aeroForces_G[2] / m;
}

void otr::calc_d_orientationVector(){
	d_orientationVector[0] = - 0.5 * (orientationVector[4] * orientationVector[1] + orientationVector[5] * orientationVector[2] + orientationVector[6] * orientationVector[3]);
	d_orientationVector[1] = 0.5 * (orientationVector[4] * orientationVector[0] - orientationVector[5] * orientationVector[3] + orientationVector[6] * orientationVector[2]);
	d_orientationVector[2] = 0.5 * (orientationVector[4] * orientationVector[3] + orientationVector[5] * orientationVector[0] - orientationVector[6] * orientationVector[1]);
	d_orientationVector[3] = 0.5 * ( - orientationVector[4] * orientationVector[2] + orientationVector[5] * orientationVector[1] + orientationVector[6] * orientationVector[0]);

	d_orientationVector[4] = aeroTorque[0] / J[0] - (((J[2] - J[1]) * orientationVector[5] * orientationVector[6]) / J[0]);
	d_orientationVector[5] = aeroTorque[1] / J[1] - (((J[0] - J[2]) * orientationVector[4] * orientationVector[6]) / J[1]);
	d_orientationVector[6] = aeroTorque[2] / J[2] - (((J[1] - J[0]) * orientationVector[4] * orientationVector[5]) / J[2]);
}


void otr::calc_d_dt(){
	calc_d_stateVectorG();
	calc_d_orientationVector();
}


///////////////////////////////////////////////////////////////////
//						Интегрирование 		 					 //
///////////////////////////////////////////////////////////////////

void otr::integ_stateVectorG(double &dt){
	for(int i = 0; i < 6; i++){
		stateVectorG[i] += d_stateVectorG[i] * dt;
	}
}

void otr::integ_orientationVector(double &dt){
	for(int i = 0; i < 7; i++){
		orientationVector[i] += d_orientationVector[i] * dt;
	}
	double norm = sqrt(orientationVector[0]*orientationVector[0] + orientationVector[1]*orientationVector[1] +
						orientationVector[2]*orientationVector[2] + orientationVector[3]*orientationVector[3]);
	for(int i = 0; i < 4; i ++){
		orientationVector[i] /= norm;
	}
}


void otr::integrate(double &dt){
	integ_stateVectorG(dt);
	integ_orientationVector(dt);
	t += dt;
}



///////////////////////////////////////////////////////////////////
//						Всё в нужном порядке	 		 		 //
///////////////////////////////////////////////////////////////////

void otr::init_step_params(){
	init_extra_params();
	init_speedAngles();
	init_angleVector();
	init_b();
	init_phi_i();
	init_forceAndTorque();
	calc_d_dt();
}


///////////////////////////////////////////////////////////////////
//						Используемые функции	 		 		 //
///////////////////////////////////////////////////////////////////

void otr::update(double &dt, int& i){
	if(i % 1000 == 0){
		/*cout << t << ' ' << stateVectorG[0] << ' ' << stateVectorG[1] << ' ' << stateVectorG[2] << ' ' << stateVectorG[3] << ' ' << stateVectorG[4] << ' ' << stateVectorG[5]
		<< ' ' << orientationVector[4] << ' ' << orientationVector[5] << ' ' << orientationVector[6] << ' ' << speedAngles[0] * 180.0 / M_PI << ' ' << speedAngles[1] * 180.0 / M_PI
		<< ' ' << angleVector[0] * 180.0 / M_PI << ' ' << angleVector[1] * 180.0 / M_PI << ' ' << angleVector[2] * 180.0 / M_PI << '\n';*/
	}
	i++;
	integrate(dt);
	init_step_params();
}

vector<double> otr::get_res(double &dt, int &i){
	bool flag;
	while(stateVectorG[1] > 0){
		update(dt, i);
	}

	/*cout << t << ' ' << stateVectorG[0] << ' ' << stateVectorG[1] << ' ' << stateVectorG[2] << ' ' << stateVectorG[3] << ' ' << stateVectorG[4] << ' ' << stateVectorG[5]
		<< ' ' << orientationVector[4] << ' ' << orientationVector[5] << ' ' << orientationVector[6] << ' ' << speedAngles[0] * 180.0 / M_PI << ' ' << speedAngles[1] * 180.0 / M_PI
		<< ' ' << angleVector[0] * 180.0 / M_PI << ' ' << angleVector[1] * 180.0 / M_PI << ' ' << angleVector[2] * 180.0 / M_PI << '\n';*/
	
	vector<double> res(3);
	for(int i = 0; i < 3; i++){
		res[i] = stateVectorG[i];
	}
	return res;
}