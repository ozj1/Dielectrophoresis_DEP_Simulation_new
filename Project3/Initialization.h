#pragma once
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
class Initialization
{
public:
	//intialization
#define pi  3.141592653589793
	double g = 9.8, Rom = 997;
	double eps0 = 8.854187817e-12, kT = 4.0453001636e-21, em = 78.5 * eps0, Di = 1;//Di is in um unit here
	//double Vpp = 40, dg = 8e-6; //for N,m,j,kg,s units, Ro density kg/m3, for N, m, j, kg, s, C units, eps0= 8.85e-12 C2/Nm2
	double acceptedtrials = 0., totaltrials = 0., acceptancerate;
	int jmax = 200000000;
	double rad, A, B, MetropolisMC, kBoltzmann = 1.38064852e-23, Temp = 293.;
	double Diameter_post = 20., post_distance = 7.; //diameter of each post, tilted distance between posts  
	double X_cyl1_offset = 10, Y_cyl1_offset = 10;
	double number_of_slides_in_comsol_output = 30;

	double post_H = 4, pH = 9.43, R = Di * (1.0e-6) / 2., eps = 78.5, e0 = 1.60217657e-19, Na = 6.022e+23, zetapost = -23.88;
	double A1w2, Aw = 3.7e-20;

	static constexpr int Particle_num = 2000;
	double paticles_cut_off_distance = 15.;

	//string input = "../../../inputfiles/di20dis7-30zSlices-albanie-2post-z=30-v=5-fine-results-mesh0.7-0.8-1MHz.txt"; //input file name
	string input = "di20dis7-30zSlices-albanie-2post-z=30-v=5-fine-results-mesh0.7-0.8-1MHz.txt"; //input file name
	string output = "out-sio2-2000particles-di20-dist7-v5-1MHz.txt";//output file name
	bool plot_dep_enrgy = false;//which type output we want
	bool Xvalues = false, Yvalues = false, Zvalues = false, DEPvalues = true;
	double Zd_cross = 6;

	string type1 = "sio2", type2 = "tio2";//particle types
	static constexpr double portion_of_the_first_particle_type = 1.;
	string SphSphElecCalType="NW_project";//SphSphElecCalType Sphere Sphere Electrostatic Calclation Type bevan or NW_project?

};

