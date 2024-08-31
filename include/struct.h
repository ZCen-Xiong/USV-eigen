#pragma once

/**
 * @FilePath     : /CSSC/include/struct.h
 * @Description  :
 * @Author       : Ye Ji 1317907830@qq.com
 * @Version      : 0.0.1
 * @LastEditors  :
 * @LastEditTime : 2024-01-02 15:02:21
 * @Copyright (c) 2023 by BIT, All Rights Reserved.
**/

#include <string>

enum ObjectType
{
	PATH = 1,
	AREA = 2,
	ENERGY = 3,
	TIME = 4
};

struct ShipPara
{
	std::string name;
	int index;
	ObjectType obj;

	int o_tar_num;

	int tar_num;
	double* tar_t;
	double* tar_x;
	double* tar_y;

	int time_num;
	double* time_commu_lb;
	double* time_commu_ub;

	double amax;
	double amin;
	double wmax;
	double wmin;

	double vcmax;
	double vcmin;
	double kcmax;
	double kcmin;
};


struct InputParas
{
	std::string name;   // 问题名称
	ObjectType mytype;  // 选择的性能指标类型

	int index;          // 问题序号

	double x1s_0;       // 开始机动调整时， x坐标 m
	double y1s_0;       // 开始机动调整时， y坐标 m
	double v1s_0;       // 开始机动调整时， 速率 m/s
	double k1s_0;       // 开始机动调整时， 航向角 rad

	double x1s_f;       // 结束机动调整时， x坐标 m 
	double y1s_f;       // 结束机动调整时， y坐标 m 
	double v1s_f;       // 结束机动调整时， 速率 m/s 
	double k1s_f;       // 结束机动调整时， 航向 rad 

	double t0;       // 开始机动调整的时间 s
	double ts;       // 开始通信的时间 s
	double te;       // 结束通信的时间 s
	double tf;       // 结束机动调整的时间 s

	double amax;     // 加速度上限 m/s^2
	double amin;     // 加速度下限 m/s^2
	double wmax;     // 航向角速度上限 rad/s
	double wmin;     // 航向角速度下限 rad/s

	double vmax;     // 机动段速度上限 m/s
	double vmin;     // 机动段速度下限 m/s
	double kmax;     // 机动段航向角上限 rad
	double kmin;     // 机动段航向角下限 rad

	double vcmax;    // 通信段速度上限 m/s
	double vcmin;    // 通信段速度下限 m/s
	double kcmax;    // 通信段航向角上限 rad
	double kcmin;    // 通信段航向角下限 rad

};