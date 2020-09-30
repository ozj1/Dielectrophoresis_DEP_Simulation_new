#pragma once
#include <string> 
#include<vector> 
#include "Initialization.h"

using namespace std;
class mutual_interactions
{
public:
	Initialization init;

	#define pi  3.141592653589793
	double eps0 = init.eps0, kT = init.kT, em = init.em, ep = init.ep, Di = init.Di;//Di is in um unit here
	double post_H = init.post_H, pH = init.pH, R = init.R, eps = init.eps, e0 = init.e0, Na = init.Na, zetasphere = init.zetasphere, zetapost = init.zetapost;
	double A1w2, Asio2 = init.Asio2, Aw = init.Aw;

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

