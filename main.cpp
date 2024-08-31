
/**
 * @FilePath     : /CSSCc:/Users/13179/Desktop/main.cpp
 * @Description  :
 * @Author       : Ye Ji 1317907830@qq.com
 * @Version      : 0.0.1
 * @LastEditors  :
 * @LastEditTime : 2024-01-02 18:10:17
 * @Copyright (c) 2023 by BIT, All Rights Reserved.
**/

#include <stdio.h>
#include <math.h>
#include "funcs.h"
#include "struct.h"


InputParas Convert2Input(ShipPara ship1) {
    InputParas input;
    input.name = ship1.name;
    input.mytype = ship1.obj;
    input.index = ship1.index;

    input.x1s_0 = ship1.tar_x[0];
    input.y1s_0 = ship1.tar_y[0];
    double vx0 = (ship1.tar_x[1] - ship1.tar_x[0]) / (ship1.tar_t[1] - ship1.tar_t[0]);
    double vy0 = (ship1.tar_y[1] - ship1.tar_y[0]) / (ship1.tar_t[1] - ship1.tar_t[0]);
    input.v1s_0 = sqrt(vx0 * vx0 + vy0 * vy0);
	printf("v_0 = %.4f\n", input.v1s_0);
    input.k1s_0 = atan2(vy0, vx0);
    input.x1s_f = ship1.tar_x[ship1.tar_num - 1];
    input.y1s_f = ship1.tar_y[ship1.tar_num - 1];
    double vxf = (ship1.tar_x[ship1.tar_num - 1] - ship1.tar_x[ship1.tar_num - 2]) / (ship1.tar_t[ship1.tar_num - 1] - ship1.tar_t[ship1.tar_num - 2]);
    double vyf = (ship1.tar_y[ship1.tar_num - 1] - ship1.tar_y[ship1.tar_num - 2]) / (ship1.tar_t[ship1.tar_num - 1] - ship1.tar_t[ship1.tar_num - 2]);
    input.v1s_f = sqrt(vxf * vxf + vyf * vyf);
    input.k1s_f = atan2(vyf, vxf);

    input.t0 = ship1.tar_t[0];
    input.ts = ship1.time_commu_lb[0];
    input.te = ship1.time_commu_ub[ship1.time_num - 1];
    input.tf = ship1.tar_t[ship1.tar_num - 1];

    input.amax = ship1.amax;
    input.amin = ship1.amin;
    input.wmax = ship1.wmax;
    input.wmin = ship1.wmin;
    input.vcmax = ship1.vcmax;
    input.vcmin = ship1.vcmin;
    input.kcmax = ship1.kcmax;
    input.kcmin = ship1.kcmin;
    // ???????????????????
	double kn = 1852.0 / 3600.0;
    input.vmax = 17 * kn;
    input.vmin = 0.5 * kn;
    double PI = 3.14159265358979323846;
    input.kmax = PI/2;
    input.kmin = -PI/2;

    return input;
}



int main()
{

    // MatEigneTest1();
    // printf("hello CSSC.\n");

    // ?????????????
    std::string _name = "ship";
    int _index = 1;
    ObjectType _obj = ENERGY;
    int _o_tar_num = 12;

	double kn = 1852.0 / 3600.0;

    int _tarnum = 5;
    double tar_t[5] = { 0.0, 525.0, 1150.0, 1675.0, 2100.0 };
    double tar_x[5] = { 0.0, 2698.5, 5397.0, 8095.5, 10794.0 };
    double tar_y[5] = { 1000.0, 1000.0, 1000.0, 1000.0, 1000.0 };

    int time_num = 2;
    double time_commu_lb[2] = { 400.0, 580.0 };
    double time_commu_ub[2] = { 490.0, 670.0 };

    double _amax = 0.04;
    double _amin = -0.04;
    double _wmax = 0.05;
    double _wmin = -0.05;

    double _vcmax = 4.0 * kn;
    double _vcmin = 2.0 * kn;
    double _kcmax = 5/6*3.14156926;
    double _kcmin = 1/6*3.14156926;


    // ????帳???
    ShipPara ship1;
    ship1.name = _name;
    ship1.index = _index;
    ship1.obj = _obj;
    ship1.o_tar_num = _o_tar_num;

    ship1.tar_num = _tarnum;
    ship1.tar_t = new double[ship1.tar_num];
    ship1.tar_x = new double[ship1.tar_num];
    ship1.tar_y = new double[ship1.tar_num];
    for (int i = 0; i < ship1.tar_num; i++) {
        ship1.tar_t[i] = tar_t[i];
        ship1.tar_x[i] = tar_x[i];
        ship1.tar_y[i] = tar_y[i];
        // ship1.tar_y[i] = -1000.0;
    }

    ship1.time_num = time_num;
    ship1.time_commu_lb = new double[ship1.time_num];
    ship1.time_commu_ub = new double[ship1.time_num];
    for (int i = 0; i < ship1.time_num; i++) {
        ship1.time_commu_lb[i] = time_commu_lb[i];
        ship1.time_commu_ub[i] = time_commu_ub[i];
    }

    ship1.amax = _amax;
    ship1.amin = _amin;
    ship1.wmax = _wmax;
    ship1.wmin = _wmin;

    ship1.vcmax = _vcmax;
    ship1.vcmin = _vcmin;
    ship1.kcmax = _kcmax;
    ship1.kcmin = _kcmin;

    // ??BUAA?????
    InputParas input1 = Convert2Input(ship1);
    // ???????
    bool flag = FeasiableSolution(input1);


    return 0;
}
