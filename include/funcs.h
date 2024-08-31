#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>
#include <functional>
#include <Eigen/Dense>
#include <random>

#include "struct.h"

/***
 * @description:
 * @param {double} theta
 * @param {double} v_comm
 * @param {double} v_max
 * @param {double} v0
 * @param {double} ts
 * @param {double} te
 * @param {double} tf
 * @param {double} acc
 * @param {double} cs
 * @param {double} ce
 * @return {*}
 */
void comm_sec(double theta, double v_comm, double v_max, double v0, double th_bias, double distance, double ts, double te, double tf, double acc, double cs[2], double ce[2]);


/***
 * @description:
 * @param {double} cart
 * @param {double} pol
 * @return {*}
 */
void cart2pol_x(double cart[2], double pol[2]);

/**
 * @brief
 * @param         {double} t0:
 * @param         {double*} y0:
 * @param         {double} h_val:
 * @param         {int} k_len:
 * @param         {double*} t:
 * @param         {double*} y:
 * @return        {*}
**/
void rk4_step(void(*odefunc)(double, double*, double*), double t0, double* y0, double h_val, int k_len, double* t, double* y);

/**
 * @brief
 * @param         {double} t0:
 * @param         {double*} y0:
 * @param         {double} h_val:
 * @param         {int} k_len:
 * @param         {double} tetol:
 * @param         {double*} t:
 * @param         {double*} y:
 * @param         {double*} te:
 * @param         {double*} ye:
 * @param         {int*} ie:
 * @return        {*}
**/

void rk4_step_events(void(*odefunc)(double, double*, double*), void(*events)(double, double*, double*, int*, int*), double t0, double* y0, double h_val, int k_len, double tetol, double* t, double* y, double* te, double* ye, int* ie);


/**
 * @brief
 * @param         {double} t:
 * @param         {double*} X:
 * @param         {double*} dX:
 * @return        {*}
**/
void dynamics_free(double t, double* X, double* dX);

/**
 * @brief
 * @param         {double} t:
 * @param         {double*} X:
 * @param         {double*} dX:
 * @return        {*}
**/
void dynamics_v_lim(double t, double* X, double* dX);

/**
 * @brief
 * @param         {double} t0:
 * @param         {double} tn:
 * @param         {double*} y0:
 * @param         {double} h:
 * @param         {int} k_len:
 * @param         {double**} t:
 * @param         {double**} y:
 * @param         {double*} te:
 * @param         {double*} ye:
 * @param         {int*} ie:
 * @return        {int} ntemp;
**/

int rk4_events(void(*odefunc)(double, double*, double*), void(*events)(double, double*, double*, int*, int*), double t0, double tn, double* y0, double h, int k_len, double* t, double** y, double* te, double* ye, int* ie);


/**
 * @brief        : �����¼�
 * @param         {double} t:
 * @param         {double*} x:
 * @param         {double*} value:
 * @param         {int*} isterminal:
 * @param         {int*} direction:
 * @return        {*}
**/
void arrival_event(double t, double* x, double* value, int* isterminal, int* direction);

/**
 * @brief        : ����ת���ĺ�����
 * @param         {double} x_arr:
 * @param         {double} xf:
 * @param         {double} x_cart:
 * @param         {int} trj_len:
 * @return        {*}
**/
void traj_plot_counter(double** x_arr, double xf, double** x_cart, int trj_len);

/**
 * @brief        : ����ת���ĺ�����
 * @param         {double**} x_arr:
 * @param         {double**} x_cart:
 * @param         {int} trj_len:
 * @return        {*}
**/
void traj_plot(double** x_arr, double** x_cart, int trj_len);

/**
 * @brief        :
 * @param         {double} th:
 * @param         {double} r:
 * @param         {double*} x:
 * @param         {double*} y:
 * @return        {*}
**/
void pol2cart(double th, double r, double* x, double* y);

/**
 * @brief        :
 * @param         {double} x:
 * @param         {double} con:
 * @return        {*}
**/
double error_continue(double x, double con[2]);

/**
 * @brief        : ������ʱ�����ſ��Ƶ���⺯��
 * @param         {double} X0:
 * @param         {double} max_TOF:
 * @param         {double} Xf:
 * @param		  {double} distance,
 * @param         {double} h:
 * @param         {int} isInverse:
 * @param         {double} e_total: ����ֵ
 * @param         {double} *TOF_opt: ����ֵ
 * @param         {double} t_arr: ����ֵ
 * @param         {double} x_f: ����ֵ
 * @param         {double} x_cart: ����ֵ
 * @return        {*}
**/
int aqua_optcon_plot(double X0[], double max_TOF, double Xf[6], double distance, double h, int isInverse, double e_total[], double* TOF_opt, double t_arr[], double x_f[], double** x_cart);

/**
 * @brief        : ������ʱ�����ſ��Ƶ���⺯��
 * @param         {double} lam0:
 * @param         {double} X0:
 * @param         {double} max_TOF:
 * @param         {double} Xf:
 * @param         {double} h:
 * @param         {double} *TOF_opt: ����ֵ
 * @param         {double} t_arr: ����ֵ
 * @param         {double} x_f: ����ֵ
 * @return        {*}
**/
int aqua_optcon_bang_lim(double lam0[], double x0[], double max_TOF, double h, double* TOF_opt, double x_f[]);

/**
 * @brief        :
 * @param         {VectorXd&} lam0:
 * @param         {VectorXd&} x0:
 * @param         {double} max_TOF:
 * @param         {VectorXd&} Xf:
 * @param         {double} h:
 * @return        {*}
**/
double ini_lam_search_bang_cart_tof(const Eigen::VectorXd& lam0, const Eigen::VectorXd& x0, double max_TOF, const Eigen::VectorXd& Xf, double h);


/**
 * @brief        :
 * @param         {int} Dimension:
 * @param         {MatrixXd&} Rmin:
 * @param         {MatrixXd&} Rmax:
 * @param         {int} Max_Gen:
 * @param         {int} Particle_Number:
 * @param         {VectorXd&} x0:
 * @param         {double} max_TOF:
 * @param         {VectorXd&} Xf:
 * @param         {double} h:
 * @return        {*}
**/
std::pair<Eigen::MatrixXd, double> CLPSO(int Dimension, const Eigen::MatrixXd& Rmin, const Eigen::MatrixXd& Rmax, int Max_Gen, int Particle_Number, const Eigen::VectorXd& x0, double max_TOF, const Eigen::VectorXd& Xf, double h);


/**
 * @brief        : whether there is a feasiable solution to the problem
 * @param         {InputParas} problem:
 * @return        {bool}: true/there is, false/none
**/
bool FeasiableSolution(InputParas problem);