#pragma once
#include <string> 
using namespace std;

class Particle
{

public:
	double eps0 = 8.854187817e-12, em = 78.5 * eps0;
	double x[2];
	double y[2];
	double z[2];
	double E[2];
	double zeta;
	double ep;
	double Asphere;
	double density;
	double CMf;
	string color;
	void set_values(double xx, double yy, double zz, double EE, int i);
    void set_particle_properties(string type);
};

