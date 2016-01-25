//
//  main.cpp
//  ToolPathTimeEstimator
//
//  Created by Celr on 15/12/10.
//  Copyright © 2015年 Celr. All rights reserved.
//
#include "ToolPathTimeEstimator.hpp"

using namespace std;
int main(int argc, const char * argv[]) {
    ToolPathTimeEstimator est = ToolPathTimeEstimator();
    char line = NULL;
    double x,y;
	int path_num;
	int line_num;
	string tt;
    size_t t = 0;
    int aa = 0;
    FILE * fd, * out, * detail, * stat;
	string name = "three_pockets_0_hollow_60.0_0.4";
	string path = "D://hollow/";
	string stat_t = ".stat";
	string out_t = ".out";
	string detail_t = ".detail";
	string dat_t = ".dat";
	string a = path + name + stat_t;
	fopen_s(&stat, (path + name + stat_t).c_str(), "w");
	fopen_s(&out, (path + name + out_t).c_str(), "w");
	fopen_s(&detail, (path + name + detail_t).c_str(), "w");
    fopen_s(&fd, (path + name + dat_t).c_str(), "r");
	aa = fscanf_s(fd, "%d", &path_num);
	for (size_t i = 0; i < path_num; i++)
	{
		aa = fscanf_s(fd, "%d", &line_num);
		for (size_t i = 0; i < line_num; i++)
		{
			
			aa = fscanf_s(fd, "%lf %lf", &x, &y);
			if (aa == 0) {
				fscanf_s(fd, "%s", &line);
			}
			else {
				if (i == 0) {
					est.addJump(x, y);
				}
				else {
					est.addBlock(x, y);
				}
			}
		}
	}

	est.prepare();


	for (int i = 20; i <= 200; i++)
	{
		std::cout << i << std::endl;
		est.nominal_speed = Vector2D(i, i);
		est.nominal_speed.direction = false;
		est.updateNominalSpeed();
		fprintf_s(stat, "%d %f\n", i, est.calculate(out));
	}

	est.nominal_speed = Vector2D(80, 80);
	est.nominal_speed.direction = false;
	est.updateNominalSpeed();
	cout << est.calculate(out) << endl;
	//est.detail(detail);
    
	est.timeline(detail, 0.01);
	//getchar();

    return 0;
}
