#pragma once
#include<vector> 
#include "MC_funcs.h"
#include "Initialization.h"

class Energy_extraction
{
public:
	Initialization init;
	double number_of_slides_in_comsol_output = init.number_of_slides_in_comsol_output;
	double Diameter_post = init.Diameter_post, post_distance = init.post_distance; //diameter of each post, tilted distance between posts  
	double dist_posts_offset = init.dist_posts_offset, X_cyl1_offset = init.X_cyl1_offset, Y_cyl1_offset = init.Y_cyl1_offset;
	struct xy
	{
		double x;
		double y;

	};
	xy coordinate_transform(double Wcell_input, double Wcell_large, double x, double y);
	double Energy_Value_Extaction(double DataNum, double DataSection, double th, double Xd, double Yd, double Zd, vector<double> X, vector<double> Y, vector<double> Z, vector<double> E, double data_lx, double data_ly, double data_lz, double d_lx);
};

