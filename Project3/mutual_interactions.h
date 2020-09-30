#pragma once
#include <string> 
#include<vector> 
using namespace std;
class mutual_interactions
{
public:
	#define pi  3.141592653589793
	double eps0 = 8.854187817e-12, kT = 4.0453001636e-21, em = 78.5 * eps0, ep = 3.8 * eps0, Di = 1;//Di is in um unit here
	double post_H = 4, pH = 9.43, R = Di * (1.0e-6) / 2., eps = 78.5, e0 = 1.60217657e-19, Na = 6.022e+23, zetasphere = -43.88, zetapost = -23.88;
	double A1w2, Asio2 = 6.6e-20, Aw = 3.7e-20;

	double ys(double zeta);

	double kappa(double pH); 

	double yy(double r, double zeta); 

	//electroestatic between sphere sphere from NW assembly paper
	//double ESt_SR_SR(double r, double zeta1, double zeta2)
	double ESt_SR_SR(double r, double zeta1, double zeta2); //zeta1 is the zeta of particle 1

	double Hr(double x, double y); 

	//double EvdwSS(double r, double Asphere1, double Asphere2)
	double EvdwSS(double r, double Asphere1, double Asphere2); //Asphere1 is the first particle type hamaker constant

	double EDL_S_FP_bevan(double h, double zetasphere); //Electric double layer between Sphere and Flat Plate Equation

	double EDL_S_FP(double h, double zetasphere); //Electric double layer between Sphere and Flat Plate Equation

	double SEIntegralSumEq(double lowBound, double upBound, int n, double z, double zetasphere);//Surface Element Integration Equation also the integral is calculated numarically 

};

