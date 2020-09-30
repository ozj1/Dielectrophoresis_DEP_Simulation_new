#pragma once
class Feild_Energy_Particle
{
public:
#define pi  3.141592653589793
	double eps0 = 8.854187817e-12, kT = 4.0453001636e-21, em = 78.5 * eps0, Di = 1;//Di is in um unit here
	double DEP(double Es, double CMf);
};

