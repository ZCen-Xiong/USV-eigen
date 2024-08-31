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
	std::string name;   // ��������
	ObjectType mytype;  // ѡ�������ָ������

	int index;          // �������

	double x1s_0;       // ��ʼ��������ʱ�� x���� m
	double y1s_0;       // ��ʼ��������ʱ�� y���� m
	double v1s_0;       // ��ʼ��������ʱ�� ���� m/s
	double k1s_0;       // ��ʼ��������ʱ�� ����� rad

	double x1s_f;       // ������������ʱ�� x���� m 
	double y1s_f;       // ������������ʱ�� y���� m 
	double v1s_f;       // ������������ʱ�� ���� m/s 
	double k1s_f;       // ������������ʱ�� ���� rad 

	double t0;       // ��ʼ����������ʱ�� s
	double ts;       // ��ʼͨ�ŵ�ʱ�� s
	double te;       // ����ͨ�ŵ�ʱ�� s
	double tf;       // ��������������ʱ�� s

	double amax;     // ���ٶ����� m/s^2
	double amin;     // ���ٶ����� m/s^2
	double wmax;     // ������ٶ����� rad/s
	double wmin;     // ������ٶ����� rad/s

	double vmax;     // �������ٶ����� m/s
	double vmin;     // �������ٶ����� m/s
	double kmax;     // �����κ�������� rad
	double kmin;     // �����κ�������� rad

	double vcmax;    // ͨ�Ŷ��ٶ����� m/s
	double vcmin;    // ͨ�Ŷ��ٶ����� m/s
	double kcmax;    // ͨ�Ŷκ�������� rad
	double kcmin;    // ͨ�Ŷκ�������� rad

};