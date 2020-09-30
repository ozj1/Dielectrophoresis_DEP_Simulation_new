#include "mutual_interactions.h"
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


double mutual_interactions::ys(double zeta) {
	return (zeta*e0) / (kT);
}
double mutual_interactions::kappa(double pH) {
	double c = pow(10., -pH)* 1.e+3;
	return sqrt((2.* pow(e0, 2.)*Na*c) / (eps* eps0* kT));
}
double mutual_interactions::yy(double r, double zeta) {
	return 4.*exp(0.5 *kappa(pH)* (r - 2.*R)) * atanh(exp(-1.*kappa(pH)*(r - 2.*R) / 2.)* tanh(ys(zeta) / 4.));
}
	//electroestatic between sphere sphere from NW assembly paper
	//double ESt_SR_SR(double r, double zeta1, double zeta2)
double mutual_interactions::ESt_SR_SR(double r, double zeta1, double zeta2) {//zeta1 is the zeta of particle 1
	r = r * 1.0e-6;
	return eps * eps0*pow((kT / e0), 2.)* yy(r, zeta1)*yy(r, zeta2)*(pow(R, 2.) / r)*log(1 + exp(-kappa(pH) * (r - 2.*R)));//modified eq for wo sphrs with differnt zeta potntial (particle spheres and post spheres) 
}
double mutual_interactions::Hr(double x, double y) {

	return y / (pow(x, 2.) + x * y + x) + y / (pow(x, 2.) + x * y + x + y) + 2.*log((pow(x, 2.) + x * y + x) / (pow(x, 2.) + x * y + x + y));
}
	//double EvdwSS(double r, double Asphere1, double Asphere2)
double mutual_interactions::EvdwSS(double r, double Asphere1, double Asphere2) {//Asphere1 is the first particle type hamaker constant
	r = r * 1.0e-6;
	A1w2 = (pow(Asphere1, 1. / 2) - pow(Aw, 1. / 2)) * (pow(Asphere2, 1. / 2) - pow(Aw, 1. / 2));
	//no problem if R and ds are not in nm scale because x and y in Hr function are ratios
	//encourted a strange problem don't write 1/12 will make zero in results, use 1./12 !!!!
	return (-1. / 12)* (A1w2*Hr((r - 2.*R) / (2.*R), 1.));
}
double mutual_interactions::EDL_S_FP_bevan(double h, double zetasphere) {//Electric double layer between Sphere and Flat Plate Equation
	//for now we consider the post zeta as -23.88
	double sai1 = tanh(ys(zetasphere) / 4.);
	double sai2 = tanh(ys(zetapost) / 4.);
	return 64 * pi *R* eps*eps0*sai1*sai2*pow((kT / e0), 2.)*exp(-kappa(pH) * ((h - post_H)* (1.0e-6) - R));//eq from SEI paper: Bhattacharjee, Elimelech
	//return r;
}
double mutual_interactions::EDL_S_FP(double h, double zetasphere) {//Electric double layer between Sphere and Flat Plate Equation
	//for now we consider the post zeta as -23.88
	double sai1 = tanh(ys(zetasphere) / 4.);
	double sai2 = tanh(ys(zetapost) / 4.);
	return 32 * eps*eps0*kappa(pH)*sai1*sai2*pow((kT / e0), 2.)*exp(-kappa(pH) * (h));//eq from SEI paper: Bhattacharjee, Elimelech
	//return r;
}
	//a function to calculate integral numarically
double mutual_interactions::SEIntegralSumEq(double lowBound, double upBound, int n, double z, double zetasphere) {//Surface Element Integration Equation also the integral is calculated numarically 
	//z = 4.01;
	n = 1000;
	double dr = (double)(upBound - lowBound) / n;
	double cumSum = 0;
	for (int i = 0; i < n; i++)
	{
		double H = (z - post_H)* (1.0e-6);//4um is the height of the posts
		//double D = H - R;
		double r = lowBound + i * dr;
		double h1 = H - R * pow((1 - pow((r / R), 2)), 0.5);
		double h2 = H + R * pow((1 - pow((r / R), 2)), 0.5);
		double funcValue = EDL_S_FP(h1, zetasphere) - EDL_S_FP(h2, zetasphere);
		double rectangleArea = funcValue * r*dr;
		cumSum += rectangleArea;

	}
	double yyy = 2 * pi*cumSum;
	return 2 * pi*cumSum;
}


