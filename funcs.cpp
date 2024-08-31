/***
 * @Author: Ye Ji 1317907830@qq.com
 * @Date: 2023-12-17 19:22:13
 * @LastEditors:
 * @LastEditTime: 2023-12-17 20:19:15
 * @FilePath: \BUAA-test\funcs.cpp
 * @Description:
 * @Copyright (c) BIT 2023, All Rights Reserved.
 */

#include "funcs.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <Eigen/Dense>
#include <iostream> 


double t_decc;
double v_max;
double v_min;
double acc_max;
double w_max;
double x_turncate;
double tf;
double v_0;
double v_f;

void comm_sec(double theta, double v_comm, double v_max, double v0,double th_bias, double distance, double ts, double te, double tf, double acc, double cs[2], double ce[2]) {
	double t_comm = te - ts;
	double d_manu_0 = v_max * (tf - t_comm) - pow((v_max - v_comm), 2) / acc - pow((v_max - v0), 2) / acc;
	double d_manu_no = fmin(d_manu_0, v0 * tf);
	double d_manu_1 = v0 * (tf - t_comm) - pow((v0 - v_comm), 2) / acc;
	double d_manu = (d_manu_no + d_manu_1) / 2;
	double d_comm_0 = v0 * tf - d_manu;
	double d_comm = d_comm_0;
	double v_comm_desire = d_comm * cos(theta) / t_comm;
	double vc;
	if (v_comm < v_comm_desire*2) {
		vc = v_comm;
	}
	else {
		vc = v_comm_desire;
	}
	double x = v0 * tf;
	double t1 = ts;
	double t_com = te - ts;
	double t2 = tf - te;
	double dis_com = vc * t_com;
	double dis_com_x = dis_com * cos(theta);
	double dis_com_y = dis_com * sin(theta);
	double non_y = -dis_com_y;
	double non_x = x - dis_com_x;
	cs[0] = non_x * t1 / (t1 + t2);
	cs[1] = non_y * t1 / (t1 + t2);
	ce[0] = x - non_x * t2 / (t1 + t2);
	ce[1] = 0 - non_y * t2 / (t1 + t2);
	// 以下为0113新添加修正项
	if (v0 < v_comm)
	{
		double x1_max = (v0 + v_comm) * ts / 2;
		double bias_x = cs[0] - x1_max * 0.8;
		cs[0] = cs[0] - bias_x;
		ce[0] = ce[0] - bias_x;
	}

	if (v0 >= v_comm)
	{
		double t_gap = (v0 - v_comm) / acc;
		double t_theory = (t_gap + ts) / 2;
		double v_max_theory = t_theory * acc + v_comm;

		if (v_max > v_max_theory)
		{
			double x1_max = (v0 + v_max_theory) * (t_theory - t_gap) / 2 + (v_comm + v_max_theory) * t_theory / 2;
			double bias_x = cs[0] - x1_max * 0.8;
			cs[0] = cs[0] - bias_x;
			ce[0] = ce[0] - bias_x;
		}
	}

	double y1_min = cs[0] * tan(th_bias);

	double bias_y = cs[1] - (y1_min *t_com / (t1 + t_com) + cs[1] * t1 / (t1 + t_com));
	cs[1] = cs[1] - bias_y;
	ce[1] = ce[1] - bias_y;
}


void cart2pol_x(double cart[2], double pol[2]) {
	double x = cart[0];
	double y = cart[1];
	double R = sqrt(pow(x, 2) + pow(y, 2));
	double al = atan2(y, x);
	pol[0] = R;
	pol[1] = al;
}


void rk4_step(void(*odefunc)(double, double*, double*), double t0, double* y0, double h_val, int k_len, double* t, double* y) {
	double* k1 = (double*)calloc(k_len, sizeof(double));
	double* k2 = (double*)calloc(k_len, sizeof(double));
	double* k3 = (double*)calloc(k_len, sizeof(double));
	double* k4 = (double*)calloc(k_len, sizeof(double));
	double* temp = (double*)calloc(k_len, sizeof(double));
	odefunc(t0, y0, k1);

	for (int i = 0; i < k_len; i++) {
		temp[i] = y0[i] + h_val * k1[i] / 2;
	}

	odefunc(t0 + h_val / 2, temp, k2);

	for (int i = 0; i < k_len; i++) {
		temp[i] = y0[i] + h_val * k2[i] / 2;
	}

	odefunc(t0 + h_val / 2, temp, k3);

	for (int i = 0; i < k_len; i++) {
		temp[i] = y0[i] + h_val * k3[i];
	}

	odefunc(t0 + h_val, temp, k4);

	*t = t0 + h_val;

	for (int i = 0; i < k_len; i++) {
		y[i] = y0[i] + h_val * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
	}

	free(k1);
	free(k2);
	free(k3);
	free(k4);
	free(temp);
}

void rk4_step_events(void(*odefunc)(double, double*, double*), void(*events)(double, double*, double*, int*, int*), double t0, double* y0, double h_val, int k_len, double tetol, double* t, double* y, double* te, double* ye, int* ie) {
	double* eventValue = (double*)calloc(k_len, sizeof(double));
	double* eventValueNext = (double*)calloc(k_len, sizeof(double));

	int terminal_useless = 1;
	int direction_useless = 0;

	events(t0, y0, eventValue, &terminal_useless, &direction_useless);
	rk4_step(odefunc, t0, y0, h_val, k_len, t, y);
	events(*t, y, eventValueNext, &terminal_useless, &direction_useless);

	if (eventValueNext[0] == 0) {
		*te = *t;
		*ye = y[0];
		*ie = 0;
	}
	else if (eventValue[0] * eventValueNext[0] < 0) {
		*te = (*t + t0) / 2;
		*ye = (y[0] + y0[0]) / 2;
		*ie = 0;
	}

	free(eventValue);
	free(eventValueNext);
}

void dynamics_free(double t, double* X, double* dX) {
	// states
	double R = X[0];
	double al = X[1];
	double v = X[2];
	double q = X[3];

	// As per the assumptions of this problem, the spacecraft has no tangential acceleration, so the thata turning rate is only reflected on q

	double dlam11 = 0; // R
	double dlam12 = 0; // alpha
	double dlam21 = 0;
	double dlam22 = 0;
	double dR = cos(q) * v;
	double dal = 0;
	double dv = 0;
	double dq = 0;

	// Populate dX
	dX[0] = dR;
	dX[1] = dal;
	dX[2] = dv;
	dX[3] = dq;
	dX[4] = dlam11;
	dX[5] = dlam12;
	dX[6] = dlam21;
	dX[7] = dlam22;
}

void dynamics_v_lim(double t, double* X, double* dX) {
	/* ???????????????
	*/
	extern double t_decc;
	extern double v_max;
	extern double v_min;
	extern double acc_max;
	extern double w_max;

	// states
	double R = X[0];
	double al = X[1];
	double v = X[2];
	double q = X[3];
	// costates
	double lam11 = X[4];  // R
	double lam12 = X[5];  // alpha
	double lam21 = X[6];  // v
	double lam22 = X[7];  // q = theta-al 

	double dlam11 = lam12 * sin(q) * v / (R * R); // R
	double dlam12 = 0;                     // alpha
	double dlam21 = -lam11 * cos(q) - lam12 * sin(q) / R;
	double dlam22 = lam11 * sin(q) * v - lam12 * cos(q) * v / R;

	double acc;
	if (lam21 > 0 && v > v_min) {
		acc = -acc_max;
		if (acc * 0.1 + v < v_min - 1e-4)
			acc = (v_min - v) * 10;
	}
	else if (lam21 < 0 && v < v_max) {
		acc = acc_max;
		if (acc * 0.1 + v > v_max + 1e-4)
			acc = (v_max - v) * 10;
	}
	else {
		acc = 0;
	}

	if (t > t_decc) {
		acc = -acc_max;
		if (acc * 0.1 + v < v_min - 1e-4)
			acc = (v_min - v) * 10;
	}

	double w;
	if (lam22 > 1) {
		w = -w_max;
	}
	else if (lam22 < -1) {
		w = w_max;
	}
	else {
		w = 0;
	}

	double dR = cos(q) * v;
	double dal = sin(q) * v / R;
	double dv = acc;
	double dq = w;

	dX[0] = dR;
	dX[1] = dal;
	dX[2] = dv;
	dX[3] = dq;
	dX[4] = dlam11;
	dX[5] = dlam12;
	dX[6] = dlam21;
	dX[7] = dlam22;
}

int rk4_events(void(*odefunc)(double, double*, double*), void(*events)(double, double*, double*, int*, int*), double t0, double tn, double* y0, double h, int k_len, double* t, double** y, double* te, double* ye, int* ie) {
	int n = (int)((tn - t0) / h);
	int ntemp = 0;

	// *t = (double *)calloc(n + 1, sizeof(double));
	// *y = (double *)calloc((n + 1) * k_len, sizeof(double));
	/*
	(*t)[0] = t0;
	for (int i = 0; i < k_len; i++) {
		(*y)[i] = y0[i];
	}
	*/
	t[0] = t0;
	for (int i = 0; i < k_len; i++) {
		y[0][i] = y0[i];
	}


	for (int i = 0; i < n; i++) {
		/*  tetol = 0.01, ??????rk4_step_events???????????????
		rk4_step_events(odefunc, events, (*t)[i], &(*y)[i * k_len], h, k_len, 0.01, &(*t)[i + 1], &(*y)[(i + 1) * k_len], te, ye, ie);
		*/
		rk4_step_events(odefunc, events, t[i], y[i], h, k_len, 0.01, &(t[i + 1]), y[i + 1], te, ye, ie);


		// ???????
		/*
		printf("y[%d] = [", i);
		for(int mm=0; mm<k_len; mm++)
		{
			printf(" %.6f,", y[i][mm]);
		}
		printf("]\n");

		printf("y[%d] = [", i+1);
		for(int mm=0; mm<k_len; mm++)
		{
			printf(" %.6f,", y[i+1][mm]);
		}
		printf("]\n");

		printf("t[%d] = %.4f.\n", i, t[i]);
		printf("t[%d] = %.4f.\n", i+1, t[i+1]);
		*/

		ntemp = i + 2; // ???????????????????????, ntemp?????

		/* ??????NULL?ж?????????? ???te, ye, ie???????洫??????????????????
		??1???????????????(??????)??
		?????2?????????*te=NULL, *ye=NULL, *ie=NULL ???????? ??????????У??????????????渳?
		??????aqua_optcon_plot???rk_events???????? ?????1????????βΣ????(*te) = -1000.0; (?????????????????????!!!)
		*/

		// if (*te != NULL) //??????????????
		if (fabs(*te + 1000.0) > 1e-6)         //(*te)?????б?????????????д???
		{
			/* ??matlab?????????????????????????????????????????t, y???, ??????????????
			   ???????·??????????????????????????????????????????????????????????????????????
			*/
			// (*t) = (double *)realloc((*t), (i + 2) * sizeof(double));
			// (*y) = (double *)realloc((*y), (i + 2) * k_len * sizeof(double));
			// break;
			/*
			???и?????????????????????????е????ntemp?? ????void rk_events????? int rk_events, ??????????е???? ntemp;
			*/
			// ????
			break;
		}
	}
	return ntemp;
}

void arrival_event(double t, double* x, double* value, int* isterminal, int* direction) {

	/* ???????????????
	*/
	extern double x_turncate;

	// x is assumed to be an array of 8 doubles
	double R = x[0];
	double al = x[1];
	*value = R * cos(al) - x_turncate;
	if (isterminal) *isterminal = 1;
	if (direction) *direction = 0;
}

void traj_plot_counter(double** x_arr, double xf, double** x_cart, int trj_len) {
	double* R_arr = (double*)calloc(trj_len, sizeof(double));
	double* al_arr = (double*)calloc(trj_len, sizeof(double));
	for (int i = 0; i < trj_len; i++) {
		R_arr[i] = x_arr[i][0];
		al_arr[i] = x_arr[i][1];
	}

	double* x = (double*)calloc(trj_len, sizeof(double));
	double* y = (double*)calloc(trj_len, sizeof(double));
	for (int i = 0; i < trj_len; i++) {
		x[i] = xf - R_arr[i] * cos(al_arr[i]);
		y[i] = R_arr[i] * sin(al_arr[i]);
	}

	for (int i = 0; i < trj_len; i++) {
		x_arr[i][3] = al_arr[i] - x_arr[i][3];
	}

	for (int i = 0; i < trj_len; i++) {
		x_cart[i][0] = x[i];
		x_cart[i][1] = y[i];
		x_cart[i][2] = x_arr[i][2];
		x_cart[i][3] = x_arr[i][3];
	}

	free(R_arr);
	free(al_arr);
	free(x);
	free(y);
}

void traj_plot(double** x_arr, double** x_cart, int trj_len) {

	double* R_arr = (double*)calloc(trj_len, sizeof(double));
	double* al_arr = (double*)calloc(trj_len, sizeof(double));

	for (int i = 0; i < trj_len; i++) {
		R_arr[i] = x_arr[i][0];
		al_arr[i] = x_arr[i][1];
	}


	double* x = (double*)calloc(trj_len, sizeof(double));
	double* y = (double*)calloc(trj_len, sizeof(double));

	for (int i = 0; i < trj_len; i++) {
		x[i] = R_arr[i] * cos(al_arr[i]);
		y[i] = R_arr[i] * sin(al_arr[i]);
	}

	for (int i = 0; i < trj_len; i++) {
		x_arr[i][3] = al_arr[i] + x_arr[i][3];
	}

	for (int i = 0; i < trj_len; i++) {
		x_cart[i][0] = x[i];
		x_cart[i][1] = y[i];
		x_cart[i][2] = x_arr[i][2];
		x_cart[i][3] = x_arr[i][3];
	}

	free(R_arr);
	free(al_arr);
	free(x);
	free(y);
}

void pol2cart(double th, double r, double* x, double* y) {
	*x = r * cos(th);
	*y = r * sin(th);
}

double error_continue(double x, double con[2]) {
	double dis_con = con[1] - con[0];
	double low = x - con[0];
	double up = x - con[1];
	double symb = 1;
	if (low < 0) {
		symb = -1;
	}
	double e = (fabs(low) + fabs(up) - dis_con) / 2 * symb;
	return e;
}

int aqua_optcon_plot(double X0[], double max_TOF, double Xf[6], double distance, double h, int isInverse, double e_total[], double* TOF_opt, double t_arr[], double x_f[], double** x_cart)
{
	// Shooting Function
	extern double tf;
	extern double v_0;

	int k_len = 8;
	double X1[8], t_temp;
	rk4_step(dynamics_free, 0.0, X0, h, k_len, &t_temp, X1);


	double tspan[2] = { t_temp, max_TOF };

	double te_useless = -1000.0;
	double ye_useless = -1000.0;
	int ie_useless = -1000;
	int n_theory = (int)((tspan[1] - tspan[0]) / h);
	double* t_array = (double*)calloc(n_theory + 1, sizeof(double));
	double** x_array = (double**)calloc(n_theory + 1, sizeof(double*));
	for (int i = 0; i < n_theory + 1; i++)
	{
		x_array[i] = (double*)calloc(k_len, sizeof(double));
	}

	int n = rk4_events(dynamics_v_lim, arrival_event, tspan[0], tspan[1], X1, h, k_len, t_array, x_array, &te_useless, &ye_useless, &ie_useless);

	//int idx = n - 1;
	//printf("y[%d,:] = [", idx);
	//for (int mm = 0; mm < k_len; mm++)
	//{
	//	printf(" %.6f,", x_array[idx][mm]);
	//}
	//printf("]\n");

	//printf("t[%d] = %.4f.\n", n - 1, t_array[n - 1]);



	if (isInverse == 1) {
		for (int i = 0; i < n; i++) {
			t_arr[i] = tf - t_array[i];
		}
		traj_plot_counter(x_array, distance, x_cart, n);
	}
	else {
		for (int i = 0; i < n; i++) {
			t_arr[i] = t_array[i];
		}
		traj_plot(x_array, x_cart, n);
	}

	//int testidx = n - 1;
	//printf("x_cart[%d,:] = [", testidx);
	//for (int mm = 0; mm < 4; mm++)
	//{
	//	printf(" %.6f,", x_cart[testidx][mm]);
	//}
	//printf("]\n");


	* TOF_opt = t_arr[n - 1];
	for (int i = 0; i < 4; i++) {
		x_f[i] = x_array[n - 1][i];
	}

	double R = x_f[0];
	double al = x_f[1];
	double v = x_f[2];
	double q = x_f[3];

	double R_con = Xf[0];
	double al_con = Xf[1];
	double v_con[2] = { Xf[2], Xf[3] };
	double th_con[2] = { Xf[4], Xf[5] };


	double xf, yf;
	pol2cart(al_con, R_con, &xf, &yf);

	double x_ach, y_ach;
	pol2cart(al, R, &x_ach, &y_ach);
	double theta = al + q;

	double e_x = xf - x_ach;
	double e_y = yf - y_ach;
	double e_v = error_continue(v, v_con);
	double e_theta = error_continue(theta, th_con);
	double pi = 3.1415926;
	double abs_e_theta = fabs(e_theta); // 计算 e_theta 的绝对值
	abs_e_theta = fmod(abs_e_theta, 2 * pi);

	e_total[0] = e_x;
	e_total[1] = e_y;
	e_total[2] = e_v;
	e_total[3] = abs_e_theta;

	// 
	for (int i = 0; i < n_theory + 1; i++) {
		free(x_array[i]);
	}
	free(t_array);
	free(x_array);
	return n;
}


int aqua_optcon_bang_lim(double lam0[], double x0[], double max_TOF, double h, double* TOF_opt, double x_f[])
{
	// Shooting Function
	extern double t_decc;
	extern double x_turncate;
	int k_len = 8;
	double X1[8], t_temp;
	double X0[8] = { x0[0], x0[1], x0[2], x0[3], lam0[0], lam0[1], lam0[2], lam0[3] };
	rk4_step(dynamics_free, 0.0, X0, h, k_len, &t_temp, X1);

	double tspan[2] = { t_temp, max_TOF };
	double te_useless = -1000.0;
	double ye_useless = -1000.0;
	int ie_useless = -1000;
	// ??????????????????????????????Ч??
	int n_theory = (int)((tspan[1] - tspan[0]) / h);
	double* t_array = (double*)calloc(n_theory + 1, sizeof(double));
	double** x_array = (double**)calloc(n_theory + 1, sizeof(double*));
	for (int i = 0; i < n_theory + 1; i++)
	{
		x_array[i] = (double*)calloc(k_len, sizeof(double));
	}

	int n = rk4_events(dynamics_v_lim, arrival_event, tspan[0], tspan[1], X1, h, k_len, t_array, x_array, &te_useless, &ye_useless, &ie_useless);

	int idx = n - 1;
	//printf("y[%d,:] = [", idx);
	//for (int mm = 0; mm < k_len; mm++)
	//{
	//	printf(" %.6f,", x_array[idx][mm]);
	//}
	//printf("]\n");

	//printf("t[%d] = %.4f.\n", n - 1, t_array[n - 1]);


	* TOF_opt = t_array[n - 1];
	for (int i = 0; i < 4; i++) {
		x_f[i] = x_array[n - 1][i];
	}

	// 
	for (int i = 0; i < n_theory + 1; i++) {
		free(x_array[i]);
	}
	free(t_array);
	free(x_array);
	return n;
}

double ini_lam_search_bang_cart_tof(const Eigen::VectorXd& lam_0, const Eigen::VectorXd& x0, double max_TOF, const Eigen::VectorXd& Xf, double h)
{
	// ??????????
	extern double t_decc;
	extern double x_turncate;

	double R_con = Xf(0);
	double al_con = Xf(1);
	double v_con[2] = { Xf(2), Xf(3) };
	double th_con[2] = { Xf(4), Xf(5) };
	double TOF_opt = 0;
	double x_f[4];
	double lam0_double[4];
	double x0_double[4];
	// 
	for (int i = 0; i < 4; ++i) {
		lam0_double[i] = lam_0(i);
		x0_double[i] = x0(i);
	}
	int n = aqua_optcon_bang_lim(lam0_double, x0_double, max_TOF, h, &TOF_opt, x_f);

	double R = x_f[0];
	double al = x_f[1];
	double v = x_f[2];
	double q = x_f[3];

	double theta = al + q;

	// desired
	double xs, ys;
	pol2cart(al_con, R_con, &xs, &ys);

	// actual
	double x_ach, y_ach;
	pol2cart(al, R, &x_ach, &y_ach);
	double pi = 3.1415926;
	double e_x = xs - x_ach;
	double e_y = ys - y_ach;
	double e_v = error_continue(v, v_con);
	double e_theta = error_continue(theta, th_con);
	double abs_e = fabs(e_theta); // 计算 e_theta 的绝对值
	abs_e = fmod(abs_e, 2 * pi);
	double e_pos = std::sqrt(std::pow(e_x, 2) + std::pow(e_y, 2));
	double J = e_pos * 100 + abs_e * 1 / pi * 180 + TOF_opt;
	return J;
}


struct Particle
{
	Eigen::Matrix<double, 1, Eigen::Dynamic> position;
	Eigen::Matrix<double, 1, Eigen::Dynamic> velocity;
	Eigen::Matrix<double, 1, Eigen::Dynamic> pBest;
	double fitness;
	double pBestFitness;
};

std::pair<Eigen::MatrixXd, double> CLPSO(int Dimension, const Eigen::MatrixXd& Rmin, const Eigen::MatrixXd& Rmax, int Max_Gen, int Particle_Number, const Eigen::VectorXd& x0, double max_TOF, const Eigen::VectorXd& Xf, double h)
{
	int ps = Particle_Number; // ??????????
	int D = Dimension;        // ??????
	int me = Max_Gen;         // ??????????

	if (Rmin.size() == 1)
	{
		Rmin.replicate(1, D);
		Rmax.replicate(1, D);
	}

	Eigen::MatrixXd mv = 0.2 * (Rmax - Rmin);
	Eigen::MatrixXd Rmin_M = Rmin.replicate(ps, 1); // ????ps*D??λ??????????????????
	Eigen::MatrixXd Rmax_M = Rmax.replicate(ps, 1);
	Eigen::MatrixXd Vmin_M = -mv.replicate(ps, 1); // ??????С?????????0.2?????????Χ
	Eigen::MatrixXd Vmax_M = -Vmin_M;

	Eigen::VectorXd w = Eigen::VectorXd::LinSpaced(me, 0.9, 0.9) - (Eigen::VectorXd::LinSpaced(me, 1, me) * (0.7 / me));
	double c1 = 0.8; // ??????????
	double c2 = 1.49; // ?????????

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dist(0.0, 1.0);

	Eigen::MatrixXd pos = Rmin_M + (Rmax_M - Rmin_M).cwiseProduct(Eigen::MatrixXd::NullaryExpr(ps, D, [&]() { return dist(gen); }));
	Eigen::MatrixXd vel = Vmin_M + (Vmax_M - Vmin_M).cwiseProduct(Eigen::MatrixXd::NullaryExpr(ps, D, [&]() { return dist(gen); }));

	std::vector<Particle> particles(ps);
	for (int i = 0; i < ps; ++i)
	{
		particles[i].position = pos.row(i);
		particles[i].velocity = vel.row(i);
		Eigen::VectorXd lam0(4);
		for (int j = 0; j < 4; j++)
		{
			lam0(j) = particles[i].position(j);
		}
		particles[i].fitness = ini_lam_search_bang_cart_tof(lam0, x0, max_TOF, Xf, h);
		particles[i].pBest = particles[i].position;
		particles[i].pBestFitness = particles[i].fitness;
	}

	double gbestval = particles[0].pBestFitness;
	int minIndex = 0;
	for (int i = 1; i < ps; ++i)
	{
		if (particles[i].pBestFitness < gbestval)
		{
			gbestval = particles[i].pBestFitness;
			minIndex = i;
		}
	}
	Eigen::MatrixXd gbest(1, D);
	gbest = particles[minIndex].pBest;


	for (int i = 1; i < me; ++i)
	{
		// printf("%d | %-9%d in PSO.\n", i, me);
		for (int k = 0; k < ps; ++k)
		{
			particles[k].velocity = w(i) * particles[k].velocity +
				c1 * (particles[k].pBest - particles[k].position).cwiseProduct(Eigen::MatrixXd::NullaryExpr(1, D, [&]() { return dist(gen); })) +
				c2 * (gbest - particles[k].position).cwiseProduct(Eigen::MatrixXd::NullaryExpr(1, D, [&]() { return dist(gen); }));

			particles[k].velocity = particles[k].velocity.cwiseMax(-mv).cwiseMin(mv);
			particles[k].position = particles[k].position + particles[k].velocity;
			particles[k].position = particles[k].position.cwiseMax(Rmin_M.row(k)).cwiseMin(Rmax_M.row(k));
			Eigen::VectorXd lam2(4);
			for (int j = 0; j < 4; j++)
			{
				lam2(j) = particles[k].position(j);
			}
			particles[k].fitness = ini_lam_search_bang_cart_tof(lam2, x0, max_TOF, Xf, h);

			if (particles[k].fitness <= particles[k].pBestFitness)
			{
				particles[k].pBest = particles[k].position;
				particles[k].pBestFitness = particles[k].fitness;
			}

			if (particles[k].pBestFitness < gbestval)
			{
				gbest = particles[k].pBest;
				gbestval = particles[k].pBestFitness;
			}
		}

	}

	return std::make_pair(gbest, gbestval);
}


bool FeasiableSolution(InputParas problem)
{
	printf("strat judgement ....\n");

	const double PI = 3.141592653589793238462643383;
	
	// some gobal variables
	v_max = problem.vmax ;                       //（节）
	v_min = problem.vcmin ;         //（节）
	acc_max = problem.amax ;        // (节/s)
	w_max = problem.wmax;            // （rad/s)
	v_0 = problem.v1s_0 ;           // （节）
	v_f = problem.v1s_f ;           // （节）
	tf = problem.tf;                 // （s）

	printf("v_max = %.4f\n", v_max);
	printf("v_min = %.4f\n", v_min);
	printf("acc_max = %.4f\n", acc_max);
	printf("w_max = %.4f\n", w_max);
	printf("v_0 = %.4f\n", v_0);
	printf("tf = %.4f\n", tf);

	// 初始位置，航向
	double cart_0[2] = { 0, 0 }; // 初始笛卡尔位置
	double cart_f[2] = { problem.x1s_f - problem.x1s_0, problem.y1s_f - problem.y1s_0 }; // 终点笛卡尔位置
	double distance = sqrt(pow((cart_f[1] - cart_0[1]), 2) + pow((cart_f[0] - cart_0[0]), 2));
	double th_bias = -atan((cart_f[1] - cart_0[1])/(cart_f[0] - cart_0[0]));
	double th_0 = problem.k1s_0 + th_bias;				// 初始笛卡尔航向               
														// 如果是倾斜航迹，加偏置，让航迹变为水平
	double th_f = - problem.k1s_f - th_bias;			// 终点笛卡尔航向的反向    

	printf("th_0 = %.4f\n", th_0);
	printf("cart_0[2] = [%.4f, %.4f].\n", cart_0[0], cart_0[1]);

	// 通讯段
	double v_comm = problem.vcmax;      // 通信段速度上限 （节）
	double th_comm = problem.kcmin + th_bias;        // 通信段角度下限
	double ts = problem.ts;
	double te = problem.te;                //
	printf("v_comm = %.4f\n", v_comm);
	printf("th_comm = %.4f\n", th_comm);
	printf("ts = %.4f\n", ts);
	printf("te = %.4f\n", te);


	// some constants 
	double h = 0.1;

	// 一些减速用的系数
	double dec_coe1 = 0.0;
	double dec_coe2 = 0.0;
	double rot_bias_1 = 0.0;
	double rot_bias_2 = 0.0;

	// 两协态矢量初值
	double lam_0[4] = { 1.5e3,-8,-750,-70 };
	double lam_f[4] = { 300,5,120,9 };

	// 调用函数进行快速确定
	double rot_1_ava, rot_1, rot_2_ava, rot_2;

	// 1. 通信段约束
	double cart_s[2];
	double cart_e[2];
	// 0113修正后的生成段
	comm_sec(th_comm, v_comm, v_max, v_0, th_bias, distance, ts, te, tf, acc_max, cart_s, cart_e);

	// 2. 求出初始点和结束点的极坐标
	double pol_0[2];
	double q_0;
	cart2pol_x(cart_0, pol_0);
	q_0 = th_0 - pol_0[1];   // 初始航向笛卡尔转极坐标
	double pol_s[2];
	cart2pol_x(cart_s, pol_s);

	// 3. 第1段机动段
	double xs_R = pol_s[0];
	double xs_al = pol_s[1];
	double xs_v[2] = { v_min, v_comm };
	double xs_th[2] = { 1.0 / 6.0 * PI + th_bias , 5.0 / 6.0 * PI + th_bias };


	// 4. 第2段机动段
	double cart_e_inv[2] = { v_0 * tf - cart_e[0], cart_e[1] };
	double ts_pseu = tf - te;
	double pol_e[2];
	cart2pol_x(cart_e_inv, pol_e);
	double xe_R = pol_e[0];
	double xe_al = pol_e[1];
	double xe_v[2] = { v_min, v_comm };
	double xe_th[2] = { -5.0 / 6.0 * PI - th_bias, -1.0 / 6.0 * PI - th_bias };


	// 5. 机动段外推的截断条件
	double x_tc = cart_s[0];
	double t4dec = ts - (v_max - v_comm) / acc_max * dec_coe1 - rot_bias_1;
	x_turncate = x_tc;
	t_decc = t4dec;

	int Dimension = 4;
	int Max_Gen = 50;
	int Particle_Number = 10;

	Eigen::MatrixXd lb(1, Dimension);
	Eigen::MatrixXd ub(1, Dimension);
	for (int i = 0; i < 4; ++i)
	{
		lb(0, i) = lam_0[i] - fabs(lam_0[i]) * 10;
		ub(0, i) = lam_0[i] + fabs(lam_0[i]) * 10;
	}
	Eigen::VectorXd X_0(4);
	X_0 << pol_0[0], pol_0[1], v_0, q_0;
	double double_useless;
	Eigen::VectorXd Xs(6);
	Xs << xs_R, xs_al, xs_v[0], xs_v[1], xs_th[0], xs_th[1];
	std::pair<Eigen::MatrixXd, double> opt_res = CLPSO(Dimension, lb, ub, Max_Gen, Particle_Number, X_0, ts, Xs, h);
	Eigen::MatrixXd gbest = opt_res.first;

	// 7. 外推得到的结果
	double X_lam_pso_s[8];
	for (int i = 0; i < 4; i++)
	{
		X_lam_pso_s[i] = X_0[i];
	}
	for (int i = 0; i < 4; i++)
	{   // gbest 类型好像不对
		X_lam_pso_s[i + 4] = gbest(i);
	}

	// 这俩好像跟前面的物理含义重复定义了，不过能跑
	double XS_total[6] = { xs_R, xs_al, xs_v[0], xs_v[1], xs_th[0], xs_th[1] };
	double XE_total[6] = { xe_R, xe_al, xe_v[0], xe_v[1], xe_th[0], xe_th[1] };

	double e_s[4];
	double TOF_opt_s;
	int n_theory = (int)((ts - 0) / h) + 1;
	int k_len = 8;
	double* t_1_arr = (double*)calloc(n_theory + 1, sizeof(double));
	double** x_cart = (double**)calloc(n_theory + 1, sizeof(double*));
	for (int i = 0; i < n_theory + 1; i++)
	{
		x_cart[i] = (double*)calloc(k_len, sizeof(double));
	}
	double x_f_s[4];
	int n = aqua_optcon_plot(X_lam_pso_s, ts, XS_total, distance, h, 0, e_s, &TOF_opt_s, t_1_arr, x_f_s, x_cart);

	printf("e_s = [%.6f, %.6f, %.6f, %.6f]\n", e_s[0], e_s[1], e_s[2], e_s[3]);
	printf("TOF_opt_s = %.6f (s)\n", TOF_opt_s);
	double t_1 = t_1_arr[n - 1];

	double x_e_tc = cart_e_inv[0];
	double t4dec_e = ts_pseu - (v_max - v_comm) / acc_max * dec_coe2 - rot_bias_2;
	x_turncate = x_e_tc;
	t_decc = t4dec_e;




	for (int i = 0; i < 4; ++i)
	{
		lb(0, i) = lam_0[i] - fabs(lam_0[i]) * 10;
		ub(0, i) = lam_0[i] + fabs(lam_0[i]) * 10;
	}
	double q_f;
	q_f = th_f - pol_0[1];   // 终点航向笛卡尔转极坐标，因为解反问题，所以和第一段初始位置一样
	Eigen::VectorXd X_f(4);
	X_f << pol_0[0], pol_0[1], v_f, q_f;


	Eigen::VectorXd Xe(6);
	Xe << xs_R, xs_al, xs_v[0], xs_v[1], xs_th[0], xs_th[1];
	std::pair<Eigen::MatrixXd, double> opt_res2 = CLPSO(Dimension, lb, ub, Max_Gen, Particle_Number, X_f, ts_pseu, Xe, h);

	Eigen::MatrixXd gbest2 = opt_res2.first;

	double X_lam_pso_e[8];
	for (int i = 0; i < 4; i++)
	{
		X_lam_pso_e[i] = X_f[i];
	}
	for (int i = 0; i < 4; i++)
	{   // gbest 数据类型的问题？
		X_lam_pso_e[i + 4] = gbest2(2);
	}

	double e_e[4];
	double TOF_opt_e;
	int n_theory2 = (int)((ts_pseu - 0) / h) + 1;
	k_len = 8;
	double* t_2_arr = (double*)calloc(n_theory2 + 1, sizeof(double));
	double** x_cart2 = (double**)calloc(n_theory2 + 1, sizeof(double*));
	for (int i = 0; i < n_theory2; i++)
	{
		x_cart2[i] = (double*)calloc(k_len, sizeof(double));
	}
	double x_f_e[4];
	int n2 = aqua_optcon_plot(X_lam_pso_e, ts_pseu, XE_total, distance, h, 1, e_e, &TOF_opt_e, t_2_arr, x_f_e, x_cart2);

	printf("e_e = [%.6f, %.6f, %.6f, %.6f]\n", e_e[0], e_e[1], e_e[2], e_e[3]);
	printf("TOF_opt_e = %.6f (s)\n", TOF_opt_e);
	double t_2 = t_2_arr[n2 - 1];
	printf("t_2  = %.6f (s)\n", t_2);
	// 转弯需要的时间
	rot_1 = fabs(e_s[3]) / w_max;
	rot_2 = fabs(e_e[3]) / w_max;
	// 剩余能用于转弯的时间
	rot_1_ava = ts - t_1;
	rot_2_ava = ts_pseu - t_2;


	// 释放内存
	for (int i = 0; i < n_theory + 1; i++)
	{
		free(x_cart[i]);
	}
	for (int i = 0; i < n_theory2 + 1; i++)
	{
		free(x_cart2[i]);
	}
	free(x_cart);
	free(x_cart2);
	free(t_1_arr);
	free(t_2_arr);

	if (rot_1_ava + rot_2_ava > (rot_2 + rot_1) && rot_1_ava > 1e-1 && rot_2_ava > 1e-1)
	{
		if (sqrt(pow(e_s[0], 2) + pow(e_s[1], 2)) < sqrt(pow(cart_s[0], 2) + pow(cart_s[1], 2)) / 9 && sqrt(pow(e_e[0], 2) + pow(e_e[1], 2)) < sqrt(pow(cart_e_inv[0], 2) + pow(cart_e_inv[1], 2)) / 9)
		{
			printf("exsit\n");
		}
		else
		{
			printf("none\n");
		}
		return true;
	}
	else
	{
		printf("none\n");
		return false;
	}

}
