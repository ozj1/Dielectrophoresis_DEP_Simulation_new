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
	double eps0 = 8.854187817e-12, kT = 4.0453001636e-21, em = 78.5 * eps0, ep = 3.8 * eps0, Di = 1;//Di is in um unit here
	double Vpp = 40, dg = 8e-6; //for N,m,j,kg,s units, Ro density kg/m3, for N, m, j, kg, s, C units, eps0= 8.85e-12 C2/Nm2
	double CMfgeneral = (ep - em) / (ep + 2 * em);//a single CMf value is needed for dep energy profile, for monte carlo calc we use the specific CMf of each particle   
	double acceptedtrials = 0., totaltrials = 0., acceptancerate;
	int jmax = 200000000;
	double rad, A, B, MetropolisMC, kBoltzmann = 1.38064852e-23, Temp = 293.;
	double Diameter_post = 20., post_distance = 30.; //diameter of each post, tilted distance between posts  
	double dist_posts_offset = 0, X_cyl1_offset = 18, Y_cyl1_offset = 18;
	double number_of_slides_in_comsol_output = 30;

	double post_H = 4, pH = 9.43, R = Di * (1.0e-6) / 2., eps = 78.5, e0 = 1.60217657e-19, Na = 6.022e+23, zetasphere = -43.88, zetapost = -23.88;
	double A1w2, Asio2 = 6.6e-20, Aw = 3.7e-20;

	static constexpr int Particle_num = 1000;
	double paticles_cut_off_distance = 15.;

	string input = "30zSlices-albanie-2post-1MHz-z=30-v=5-normal-results-mesh0.6-0.8.txt"; //input file name
	string output = "Output/out sio2 2500 particles di 20 dist 30 v5 1MHzmnmn.txt";//output file name
	bool plot_dep_enrgy = false;//which type output we want
	bool Xvalues = false, Yvalues = false, Zvalues = false, DEPvalues = true;
	double Zd_cross = 6;

	string type1 = "sio2", type2 = "tio2";//particle types
	static constexpr double portion_of_the_first_particle_type = 0.5;
};

