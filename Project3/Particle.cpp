#include "Particle.h"
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

	
void Particle::set_values(double xx, double yy, double zz, double EE, int i) {
	x[i] = xx;
	y[i] = yy;
	z[i] = zz;
	E[i] = EE;
}
	
void Particle::set_particle_properties(string type) {
	if (type == "pmma") {//it means
		color = "[1, 0.5, 0.5]";//red color 
		zeta = -13.64;//mV
		ep = 2.6 * eps0;
		Asphere = 7.1e-20;//J
		density = 1.19e6;//gr/m3
		CMf = (ep - em) / (ep + 2 * em);
	}
	else if (type == "sio2") {
		//color = "[224, 255, 255]";//kinda blue
		color = "[0.6, 0.6, 1]";// light blue
		zeta = -43.88;
		ep = 3.8 * eps0;
		Asphere = 6.6e-20;
		density = 2.2e6;//gr/m3
		CMf = (ep - em) / (ep + 2 * em);
	}
	else if (type == "tio2") {
		color = "[0, 255, 0]";//green
		zeta = -39.99;
		ep = 80.* eps0;
		Asphere = 43.e-20;
		density = 3.9e6;
		CMf = (ep - em) / (ep + 2 * em);
	}

}


