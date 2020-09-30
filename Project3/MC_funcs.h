#pragma once
#include<vector> 
#include <string> 
using namespace std;

class MC_funcs
{
public:
	typedef struct {
		vector<double> xvalue;
		vector<double> yvalue;
		vector<double> zvalue;
		vector<double> Enorm;
		//int totalnum; //this is only added to save the total number of data for being used out of function
	} triple;
	triple getdata(string input);
	//i made this function because rand is useful for int variables, RAND_MAX is known in library of <stdlib.h>
	double fRand(double fMin, double fMax);
	// exactly here we can check the risk of overlap for rectangles with PBC situation, then in overlapchecker we need to do some modifications to reach our goal in PBC
	double calcDistance(double x1, double x2, double y1, double y2, double z1, double z2, double lx, double ly);
};

