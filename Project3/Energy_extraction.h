#pragma once
#include<vector> 
#include "MC_funcs.h"

class Energy_extraction
{
public:
	double number_of_slides_in_comsol_output = 30;
	double Diameter_post = 20., post_distance = 30.; //diameter of each post, tilted distance between posts  
	double dist_posts_offset = 0, X_cyl1_offset = 18, Y_cyl1_offset = 18;
	struct xy
	{
		double x;
		double y;

	};
	xy coordinate_transform(double Wcell_input, double Wcell_large, double x, double y);
	double Energy_Value_Extaction(double DataNum, double DataSection, double th, double Xd, double Yd, double Zd, vector<double> X, vector<double> Y, vector<double> Z, vector<double> E, double data_lx, double data_ly, double data_lz, double d_lx);
};

