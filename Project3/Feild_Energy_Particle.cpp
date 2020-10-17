#include "Feild_Energy_Particle.h"
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

double Feild_Energy_Particle::DEP(double Es, double CMf)
{
	double landa;// , EE0;
	//EE0 = pow(8., -0.5)*Vpp / dg;

	landa = pi * em*pow((Di*1e-6 / 2.), 3.)*pow(CMf, 2.) / kT;
	return -2. * kT*landa*pow(Es, 2.) / CMf;
	//return 2. * pi * em*CMf *pow((Di*1e-6 / 2.), 3.)*pow(Es, 2.);
}
double Feild_Energy_Particle::buoyant_force(double z, double Rop)
{
	z = z * 1.e-6;// *10 ^ -6 to account for z being in um unit
	Di = Di * 1.e-6;
	return (4 / 3.)*pi*pow((Di / 2.), 3.)*(Rop - Rom)*g*z;
}
