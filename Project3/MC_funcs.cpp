#include "MC_funcs.h"
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream> 
#include <string> 
#include <stdlib.h> 
#include <time.h> 
#include <math.h> 
#include<vector> 
#include <iomanip>      
using namespace std;
#include <algorithm> 


//struct wave_parameters linear_dispersion () {

    MC_funcs::triple MC_funcs::getdata(string input) {
	MC_funcs::triple fromprevious;
	vector<double> values;
	double in = 0.0;

	ifstream infile;
	infile.open(input);//remember to save a text file in exact location that "out" file exist  
	while (infile >> in) {
		values.push_back(in);
	}
	int jj;// so imp you should update nn each time (which is the number points in each rectangle) NN is the number of particles
	for (jj = 0; jj < int(values.size()); jj++) {
		//cout << values[j] << endl;
		if (jj % 4 == 0) {
			fromprevious.xvalue.push_back(values[jj]);
			//fromprevious.xvalue[jj / 4] = values[jj];
		}
		if (jj % 4 == 1) {
			fromprevious.yvalue.push_back(values[jj]);
			//fromprevious.yvalue[jj / 4] = values[jj];
		}
		if (jj % 4 == 2) {
			fromprevious.zvalue.push_back(values[jj]);
			//fromprevious.zvalue[jj / 4] = values[jj];// it rotates periodiacally from x-y-z
		}
		if (jj % 4 == 3) {
			fromprevious.Enorm.push_back(values[jj]);
			//fromprevious.Enorm[jj / 4] = values[jj];// E:V/m value is savd here
		}
	}
	//fromprevious.totalnum = jj;

	return fromprevious;
}

	//i made this function because rand is useful for int variables, RAND_MAX is known in library of <stdlib.h>
double MC_funcs::fRand(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}
	// exactly here we can check the risk of overlap for rectangles with PBC situation, then in overlapchecker we need to do some modifications to reach our goal in PBC
double MC_funcs::calcDistance(double x1, double x2, double y1, double y2, double z1, double z2, double lx, double ly)
{
	//nwe added: previous problem was i didn't add abs
	// no one has overlap problem
	if (abs(x1 - x2) <= (lx / 2) && abs(y1 - y2) <= (ly / 2)) {
		return sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2) + pow((z1 - z2), 2));
	}//x has overlap
	else if (abs(x1 - x2) >= (lx / 2) && abs(y1 - y2) <= (ly / 2)) {
		return sqrt(pow((lx - abs(x1 - x2)), 2) + pow((y1 - y2), 2) + pow((z1 - z2), 2));
	}//x,y has overlap
	else if (abs(x1 - x2) >= (lx / 2) && abs(y1 - y2) >= (ly / 2)) {
		return sqrt(pow((lx - abs(x1 - x2)), 2) + pow((ly - abs(y1 - y2)), 2) + pow((z1 - z2), 2));
	}//y,z has overlap
	else if (abs(x1 - x2) <= (lx / 2) && abs(y1 - y2) >= (ly / 2)) {
		return sqrt(pow((x1 - x2), 2) + pow((ly - abs(y1 - y2)), 2) + pow((z1 - z2), 2));
	}
}


