//Monte Carlo code for DEP posts work
//in here we don't use sei integration method to calculte sphere plate electro energy, we use the equation provided in Bevan's paper
//In here we account for coordination transform as comsol simulation cannot handle more than 2 cylinder posts 
//this code has been adjusted for comsol slice data  
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream> // for getting ofstream out("out.txt"); as a cout as a trick for avoiding high amount for dimensions because randNumx[plb][plb] will give us stackoverflow error
#include <string> // for getting ofstream out("out.txt");
#include <stdlib.h> // for using rand_s
// in order to using cin, cout, endl, ...
#include <time.h> // to use srand, time
#include <math.h> // because of POW sqrt, abs, ...
// in order to using cin, cout, endl, ...
#include<vector> //this code use vextor library which is said to be the most reliable way for pulling doubles from .txt
#include <iomanip>      // std::setw
using namespace std;
#include <algorithm>    // std::min_element, std::max_element

#include "Initialization.h"

//adding classes
#include "MC_funcs.h" 
#include "Particle.h" 
#include "Feild_Energy_Particle.h"
#include "Energy_extraction.h"
#include "mutual_interactions.h"

Initialization init;
	//intialization
#define pi  3.141592653589793
double eps0 = init.eps0, kT = init.kT, em = init.em, Di = init.Di;//Di is in um unit here
//double Vpp = init.Vpp, dg = init.dg; //for N,m,j,kg,s units, Ro density kg/m3, for N, m, j, kg, s, C units, eps0= 8.85e-12 C2/Nm2
double acceptedtrials = init.acceptedtrials, totaltrials = init.totaltrials, acceptancerate;
int jmax = init.jmax;
double rad, A, B, MetropolisMC, kBoltzmann = init.kBoltzmann, Temp = init.Temp;
double Diameter_post = init.Diameter_post, post_distance = init.post_distance; //diameter of each post, tilted distance between posts  
double X_cyl1_offset = init.X_cyl1_offset, Y_cyl1_offset = init.Y_cyl1_offset;
double number_of_slides_in_comsol_output = init.number_of_slides_in_comsol_output;

double post_H = init.post_H, pH = init.pH, R = init.R, eps = init.eps, e0 = init.e0, Na = init.Na, zetapost = init.zetapost;
double A1w2, Aw = init.Aw;

constexpr int Particle_num = init.Particle_num;
double paticles_cut_off_distance = init.paticles_cut_off_distance;

string input = init.input; //input file name
string output = init.output;//output file name
bool plot_dep_enrgy = init.plot_dep_enrgy;//which type output we want
bool Xvalues = init.Xvalues, Yvalues = init.Yvalues, Zvalues = init.Zvalues, DEPvalues = init.DEPvalues;
double Zd_cross = init.Zd_cross;

string type1 = init.type1, type2 = init.type2;//particle types
constexpr double portion_of_the_first_particle_type = init.portion_of_the_first_particle_type;
string SphSphElecCalType=init.SphSphElecCalType;//SphSphElecCalType Sphere Sphere Electrostatic Calclation Type bevan or NW_project?

//const int numvalue = 565820;//as we know our data has this much row in it 565820


int main() {

	//Based on statistical thermodynamics course (central limit theorem) ensemble average over time and position is the same, check it, so if you have random distribution over time steps it would be enough 
		//Because of the above reason we generate a large set of numbers in normal distribution using Box-Miller approach

	int i = 0, o = 0, zero = 0;
	double data_lx, data_ly, data_lz;
	//if you have lj then you should decrease h to 0.00000005 to avoid overlap and shooting
	//vector<double>  X, Y, Z, E;
	vector<double> X, Y, Z, E;
	//creating an iterator for the vector
	//vector<int>::iterator it;

	//new point understood, if our data doesn't have enough resolution then we gt the same value of Ed for all Z values which can lead to inaccurate simulation results 
	ofstream file(output);
	//streambuf *coutbuf = std::cout.rdbuf();
	cout.rdbuf(file.rdbuf());


	//making the initial position of particles
	srand(time(NULL));

	//int totnum, xvalnum;
	//totnum = xyz.totalnum;
	//xvalnum = int(totnum/4)+10;
	//N = xvalnum;
	Energy_extraction energy_extraction;
	mutual_interactions mutual_interactions;
	Feild_Energy_Particle Feild_Energy_Particle;
	MC_funcs MC_funcs;
	//for (it = X.begin(); it != X.end(); ++it)	{	
		//X.emplace(it + 5, xyz.xvalue[it]);
		//X[i] = xyz.xvalue[i];
		//Y[i] = xyz.yvalue[i];
		//Z[i] = xyz.zvalue[i];
		//E[i] = xyz.Enorm[i];
	//}

	MC_funcs::triple xyz;
	xyz=MC_funcs.getdata(input);
	//triple xyz;
	//xyz = getdata();


	int DataNum = int(xyz.xvalue.size());

	data_lx = xyz.xvalue[DataNum - 1], data_ly = data_lx, data_lz = *max_element(xyz.zvalue.begin(), xyz.zvalue.end());//* to get the value of max not index;
	int DataSection = int(DataNum / data_lx);//140062/70.71=1981  |3781=567224/150, 565806 are the total # of data and we know x coordinate is sorted, by doing this we can make the code faster


	double Diameter = Diameter_post;
	double Distance = post_distance;
	double bb = (Diameter + Distance) / 1.414;
	double d_lx = (Diameter + Distance) * 2 * pow(2, 0.5);//d_lx = dsired - lx meaning the large sim box;
	double d_ly = d_lx, d_lz = data_lz;
	for (int i = 0; i < DataNum; i++) {
		X.push_back(xyz.xvalue[i]);
		Y.push_back(xyz.yvalue[i]);
		Z.push_back(xyz.zvalue[i]);
		E.push_back(xyz.Enorm[i]);
	}


	double Xd = 30, Yd = 30, Zd = 5, z0 = 0., th = 3.;//th=threshold, Xd=Xdesird, Yd=Ydesired, Zd=Zdesired

	double Ed = 0.0;//Ed=Ederivd
	Ed = energy_extraction.Energy_Value_Extaction(DataNum, DataSection, th, Xd, Yd, Zd, X, Y, Z, E, data_lx, data_ly, data_lz, d_lx);

	constexpr int Pnum = Particle_num;
	//constexpr double particle_ratio = portion_of_the_first_particle_type;
	constexpr double particle_ratio = portion_of_the_first_particle_type;

	//creating paticle objects
	int const Pnumtype1 = Pnum * particle_ratio;// constexpr for const expressions 
	//Create a course of students
	//Particle Ps[Pnumtype1];//Ps: particles
	Particle Ps[Pnum];//Ps: particles

	//assigning poperties to two different prticle types making it suitable for binary assembly
	for (i = 0; i < Pnum; i++) {
		if (i < Pnumtype1) {
			Ps[i].set_particle_properties(type1);
		}
		else {
			Ps[i].set_particle_properties(type2);
		}
	}


	int k = 0, ke = 0;
	//double Px[Pnum][2], Py[Pnum][2], Pz[Pnum][2], PEnergy[Pnum][2], A, B, A1, A2, A3, A4, A5, A6;//Di is in 1um
	double A, B, A1, A2, A3, A4, A5, A6, A7, A8;//Di is in 1um
	bool inPost=false;// if inpost true then it means the particle in inside the post which shold not happen a the posts are solid, so we reduce i-- and try find a new random position for the particle
	bool OnPost = false;//OnPost = true if the particle is on top of a post 
	int OverlapChance, ff = 0;
	double G;
	if (plot_dep_enrgy == false) {
		for (i = 0; i < Pnum; i++) {
			Ps[i].x[0] = (double)MC_funcs.fRand(0., d_lx);
			Ps[i].y[0] = (double)MC_funcs.fRand(0., d_ly);
			//debugged: in direction of Z box will start from 0 upto lz not from -lz/2 upto lz/2
			//new debugged: based on exp(-kappa * (z - (Di / 2.))) z cannot be lower than radius of particle 'cause it will pass the beneath wall which is not real
			Ps[i].z[0] = (double)MC_funcs.fRand(z0 + 0.0001, d_lz);

			if (Ps[i].z[0] <= (post_H + Di / 2.)) {
				inPost = MC_funcs.DistfromPosts(Ps[i].x[0], Ps[i].y[0], Ps[i].z[0], d_lx, d_ly, bb, Diameter);
			}
	        if (inPost==false) {
				for (k = (i - 1); k >= 0; k--) {

					int OverlapChance;
					A = MC_funcs.calcDistance(Ps[i].x[0], Ps[k].x[0], Ps[i].y[0], Ps[k].y[0], Ps[i].z[0], Ps[k].z[0], d_lx, d_ly);
					//B = 1.05*Di; we make it dimensionless
					B = 1.05;
					if (A > B) {
						OverlapChance = false;
					}
					else {
						OverlapChance = true;
					}

					if (OverlapChance == true) {
						i--;
						break;
					}
				}
			}
			else if (inPost){
				inPost = false;
				i--;
			}
		}

		for (k = 0; k < Pnum; k++) {
			G = (double)zero;
			//first we need to extract enrgy value of each point from our text data
			Xd = Ps[k].x[0], Yd = Ps[k].y[0], Zd = Ps[k].z[0];
			Ed = energy_extraction.Energy_Value_Extaction(DataNum, DataSection, th, Xd, Yd, Zd, X, Y, Z, E, data_lx, data_ly, data_lz, d_lx);
			//now we add Electrostatic among particles
			for (ke = 0; ke < Pnum; ke++) {
				if (ke == k) {
					ke++;
				}
				A = MC_funcs.calcDistance(Ps[k].x[0], Ps[ke].x[0], Ps[k].y[0], Ps[ke].y[0], Ps[k].z[0], Ps[ke].z[0], d_lx, d_ly);
				B = paticles_cut_off_distance;
				// new added: we need to add A > 340 because if there's a overlab the energy of system will increase to infinity						
				if (A <= B) {
					//G += EvdwSS(A) + ESt_SR_SR(A) + ESt_SR_ER(A) + ESt_ER_ER(A);
					//G[i][1] += EvdwSS(A) + ESt_SR_SR(A);								
					//G += mutual_interactions.EvdwSS(A, Ps[k].Asphere, Ps[ke].Asphere) + mutual_interactions.ESt_SR_SR(A, Ps[k].zeta, Ps[ke].zeta);
					//G = mutual_interactions.EvdwSS(A, Ps[k].Asphere, Ps[ke].Asphere);
					//G =mutual_interactions.EDL_SR_SR_bevan(A, Ps[k].zeta, Ps[ke].zeta);
					//G = mutual_interactions.ESt_SR_SR(A, Ps[k].zeta, Ps[ke].zeta);
					if (SphSphElecCalType == "bevan") {
						G += mutual_interactions.EvdwSS(A, Ps[k].Asphere, Ps[ke].Asphere) + mutual_interactions.EDL_SR_SR_bevan(A, Ps[k].zeta, Ps[ke].zeta);
					}
					else if (SphSphElecCalType == "NW_project") {
						G += mutual_interactions.EvdwSS(A, Ps[k].Asphere, Ps[ke].Asphere) + mutual_interactions.ESt_SR_SR(A, Ps[k].zeta, Ps[ke].zeta);
					}
				}
			}
			OnPost = MC_funcs.DistfromPosts(Ps[k].x[0], Ps[k].y[0], Ps[k].z[0], d_lx, d_ly, bb, Diameter);
			if (OnPost) {
				double mbmb= mutual_interactions.EDL_S_FP_bevan(Zd, Ps[k].zeta);
				if (Ps[k].z[0] >= (post_H + Di / 2.)) {
					G += mutual_interactions.EDL_S_FP_bevan(Zd, Ps[k].zeta);
				}
			}
			else if (OnPost == false) {
				double mbmb = mutual_interactions.EDL_S_FP_bevan(Zd + post_H, Ps[k].zeta);
				G += mutual_interactions.EDL_S_FP_bevan(Zd + post_H, Ps[k].zeta);
			}
			double mdmd = Feild_Energy_Particle.buoyant_force(Zd, Ps[k].density);
			Ps[k].E[0] = Feild_Energy_Particle.DEP(Ed, Ps[k].CMf) + Feild_Energy_Particle.buoyant_force(Zd, Ps[k].density) + G;
		}

		int j = 0;
		// NU is a counter to track usage of norm[] array and if it end, so we renew it
		for (j = 0; j < jmax; j++) {
			totaltrials++;
			int min = 0, max = Pnum - 1;

			int randi = rand() % (max - min + 1) + min;

			double radnum, randdeltax, randdeltay, randdeltaz;
			radnum = (double)MC_funcs.fRand(0., 1.);
			if (0. <= radnum && radnum < (1. / 3)) {

				randdeltax = (double)MC_funcs.fRand(-d_lx / 2.5, d_lx / 2.5);
				Ps[randi].x[1] = Ps[randi].x[0] + randdeltax;
				Ps[randi].y[1] = Ps[randi].y[0];
				Ps[randi].z[1] = Ps[randi].z[0];

			}
			else if ((1. / 3) <= radnum && radnum < (2. / 3)) {

				randdeltay = (double)MC_funcs.fRand(-d_ly / 2.5, d_ly / 2.5);
				Ps[randi].x[1] = Ps[randi].x[0];
				Ps[randi].y[1] = Ps[randi].y[0] + randdeltay;
				Ps[randi].z[1] = Ps[randi].z[0];

			}
			else {
				randdeltaz = (double)MC_funcs.fRand(-(d_lz - z0) / 5., (d_lz - z0) / 5.);
				Ps[randi].x[1] = Ps[randi].x[0];
				Ps[randi].y[1] = Ps[randi].y[0];
				Ps[randi].z[1] = Ps[randi].z[0] + randdeltaz;
			}

			//pbc
			//pbc for x direction
			if (Ps[randi].x[1] > d_lx) {

				Ps[randi].x[1] = Ps[randi].x[1] - d_lx;
			}
			else if (Ps[randi].x[1] < 0.) {

				Ps[randi].x[1] = Ps[randi].x[1] + d_lx;
			}
			//pbc for y direction
			if (Ps[randi].y[1] > d_ly) {

				Ps[randi].y[1] = Ps[randi].y[1] - d_ly;
			}
			else if (Ps[randi].y[1] < 0.) {

				Ps[randi].y[1] = Ps[randi].y[1] + d_ly;
			}

			//pbc for z direction
			if (Ps[randi].z[1] > d_lz) {

				//Pz[randi][1] = lz;
				Ps[randi].z[1] = Ps[randi].z[1] - (d_lz - z0);
			}
			else if (Ps[randi].z[1] < z0) {
				//this should not happen as we have always reulsion with the below surface 
				//Pz[randi][1] = z0+0.001;
				Ps[randi].z[1] = Ps[randi].z[1] + (d_lz - z0);

			}
			if (Ps[randi].z[1] <= (post_H + Di / 2.)) {
				inPost = MC_funcs.DistfromPosts(Ps[randi].x[1], Ps[randi].y[1], Ps[randi].z[1], d_lx, d_ly, bb, Diameter);
			}
			if (inPost == false) {
				o = 0;
				for (i = 0; i < Pnum; i++) {

					ff = 0;
					if (i != randi) {
						A = MC_funcs.calcDistance(Ps[randi].x[1], Ps[i].x[0], Ps[randi].y[1], Ps[i].y[0], Ps[randi].z[1], Ps[i].z[0], d_lx, d_ly);
						//B = 1.05*Di; we make it dimensionless
						B = 1.005;
						if (A > B) {
							OverlapChance = false;
						}
						else {
							OverlapChance = true;
						}

						if (OverlapChance == true) {
							j--;
							ff = 0;
							break;
						}
						o++;

						if (o == (Pnum - 1)) {
							G = (double)zero;
							//Px[randi][1]= (d_lx / 1.)-3, Py[randi][1]= (d_ly / 1.)-3, Pz[randi][1]=4.01;
							Xd = Ps[randi].x[1], Yd = Ps[randi].y[1], Zd = Ps[randi].z[1];
							//now a loop to find PEnergy value by searching through the Edelta2 values 
							Ed = energy_extraction.Energy_Value_Extaction(DataNum, DataSection, th, Xd, Yd, Zd, X, Y, Z, E, data_lx, data_ly, data_lz, d_lx);
							for (ke = 0; ke < Pnum; ke++) {
								if (ke != randi) {
									A = MC_funcs.calcDistance(Ps[randi].x[1], Ps[ke].x[1], Ps[randi].y[1], Ps[ke].y[1], Ps[randi].z[1], Ps[ke].z[1], d_lx, d_ly);
									B = paticles_cut_off_distance;
									// new added: we need to add A > 340 because if there's a overlab the energy of system will increase to infinity						
									if (A <= B) {
										//G += EvdwSS(A) + ESt_SR_SR(A) + ESt_SR_ER(A) + ESt_ER_ER(A);
										//G[i][1] += EvdwSS(A) + ESt_SR_SR(A);								
										//G += mutual_interactions.ESt_SR_SR(A, Ps[randi].zeta, Ps[ke].zeta) + mutual_interactions.EvdwSS(A, Ps[randi].Asphere, Ps[ke].Asphere);
										if (SphSphElecCalType == "bevan") {
											G += mutual_interactions.EDL_SR_SR_bevan(A, Ps[randi].zeta, Ps[ke].zeta) + mutual_interactions.EvdwSS(A, Ps[randi].Asphere, Ps[ke].Asphere);
										}
										else if (SphSphElecCalType == "NW_project") {
											G += mutual_interactions.ESt_SR_SR(A, Ps[randi].zeta, Ps[ke].zeta) + mutual_interactions.EvdwSS(A, Ps[randi].Asphere, Ps[ke].Asphere);
										}

									}
								}
							}
							OnPost = MC_funcs.DistfromPosts(Ps[randi].x[1], Ps[randi].y[1], Ps[randi].z[1], d_lx, d_ly, bb, Diameter);
							if (OnPost) {
								//double Gu= mutual_interactions.SEIntegralSumEq(0, R, 500, Zd);
								//double Guu = mutual_interactions.EDL_S_FP_bevan(Zd);
								if (Ps[randi].z[1] >= (post_H + Di / 2.)) {

									G += mutual_interactions.EDL_S_FP_bevan(Zd, Ps[randi].zeta);
								}
							}
							else if (OnPost == false) {
								G += mutual_interactions.EDL_S_FP_bevan(Zd+ post_H, Ps[randi].zeta);
							}

							Ps[randi].E[1] = Feild_Energy_Particle.DEP(Ed, Ps[randi].CMf)+ Feild_Energy_Particle.buoyant_force(Zd, Ps[randi].density) + G;

							if (Ps[randi].E[0] >= Ps[randi].E[1]) {
								// false == 0 and true = !false
								MetropolisMC = true;
							}

							else {
								double P = exp(-(Ps[randi].E[1] - Ps[randi].E[0]) / (kBoltzmann*Temp));
								//double P = exp(-(TotalEnergy[1] - TotalEnergy[0]) / (kBoltzmann*Temp));
								rad = (double)MC_funcs.fRand(0., 1.);
								//Yes, the random number r should be less than or equal to p = exp(-Delta E/kT). This is right.
								if (P >= rad) {
									MetropolisMC = true;
								}
								else {

									MetropolisMC = false;
								}

							}


							//checking having overlap with all previous made rectangles

							if (MetropolisMC == true)
							{
								acceptedtrials++;

								//if (Px[randi][0] == Px[randi][1] && Py[randi][0] == Py[randi][1] && Pz[randi][0] == Pz[randi][1] && radnum < (2. / 3)) {

								//	cout << "";
								//}

								Ps[randi].x[0] = Ps[randi].x[1];
								Ps[randi].y[0] = Ps[randi].y[1];
								Ps[randi].z[0] = Ps[randi].z[1];
								Ps[randi].E[0] = Ps[randi].E[1];




								//cout the new positions
								acceptancerate = acceptedtrials / totaltrials;
								//cout << "\n\n " << "acceptance rate: " << acceptancerate << "\n\n";
								int pr;
								//to get monte carlo steps per each sphere
								pr = fmod(j, (Pnum));
								if (pr == 0) {
									cout << "j=" << j << "\n\n";
									for (i = 0; i < Pnum; i++) {
										//each particle has a specific color 
										cout << "{RGBColor" << Ps[i].color << ", Opacity[0.8], Sphere[{" << Ps[i].x[0] << ", " << Ps[i].y[0] << ", " << Ps[i].z[0] << "}, " << (Di / 2.) << "]},";
									}
								}
								//if (pr == 0) {
								//	gh++;
								//	cout << (Pnum) << "\n " << gh << "\n";
								//	for (i = 0; i < Pnum; i++) {

								//		cout << " P" << (i+1) << " " << Px[i][0] << "    " << Py[i][0] << "    " << Pz[i][0] << "\n";

								//	}

								//}
							}
							else {
								j--;
							}
						}
					}
				}//here
			}
			else if (inPost) {
			inPost = false;
			j--;
			}

		}
	}
	if (plot_dep_enrgy == true) {
		double CMfgeneral = Ps[0].CMf;//a single CMf value is needed for dep energy profile, for monte carlo calc we use the specific CMf of each particle   
		//Ps[0].CMf i th CMf of particle type 1 whatever it is sio2, pmma or tio2
		//this part is just for plotting dep energy
		int j = 0;
		//out.open(generatefilename("Xvalues"));
		double deltax = 2.;
		double deltay = 2.;
		Zd = Zd_cross;//the cross of Z that we want to see the results in

		if (Xvalues == true) {
			for (i = 0; (i*deltax) <= d_lx; i++) {
				for (j = 0; (j*deltay) <= d_ly; j++) {
					cout << i * deltax << "\n";
				}
			}
		}
		//out.close();

		else if (Yvalues == true) {
			//out.open(generatefilename("Yvalues"));
			for (i = 0; (i*deltax) <= d_lx; i++) {
				for (j = 0; (j*deltay) <= d_ly; j++) {
					cout << j * deltay << "\n";
				}
			}
		}
		//out.close();
		if (Zvalues == true) {

			//out.open(generatefilename("Zvalues"));
			for (i = 0; (i*deltax) <= d_lx; i++) {
				for (k = 0; (k*deltay) <= d_ly; k++) {
					cout << Zd << "\n";
				}
			}
		}

		//out.close();
		if (DEPvalues == true) {

			//out.open(generatefilename("DEPvalues"));
			for (i = 0; (i*deltax) <= d_lx; i++) {
				for (j = 0; (j*deltay) <= d_ly; j++) {
					Xd = i * deltax, Yd = j * deltay;
					Ed = energy_extraction.Energy_Value_Extaction(DataNum, DataSection, th, Xd, Yd, Zd, X, Y, Z, E, data_lx, data_ly, data_lz, d_lx);
					cout << Feild_Energy_Particle.DEP(Ed, CMfgeneral) << "\n";
				}
			}
		}
	}
	return 0;
}