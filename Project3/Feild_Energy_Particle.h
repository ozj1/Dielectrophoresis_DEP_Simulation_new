#pragma once
#include "Initialization.h"

class Feild_Energy_Particle
{
public:
	Initialization init;
#define pi  3.141592653589793
	double eps0 = init.eps0, kT = init.kT, em = init.em, Di = init.Di, g = init.g, Rom = init.Rom;//Di is in um unit here
	double DEP(double Es, double CMf);
	double buoyant_force(double z, double Rop);
};

