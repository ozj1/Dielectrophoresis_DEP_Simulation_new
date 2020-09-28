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

//intialization
#define pi 3.141592653589793
double eps0 = 8.854187817e-12, kT = 4.0453001636e-21, em = 78.5 * eps0, ep = 3.8 * eps0, Di = 1;//Di is in um unit here
double Vpp = 40, dg = 8e-6; //for N,m,j,kg,s units, Ro density kg/m3, for N, m, j, kg, s, C units, eps0= 8.85e-12 C2/Nm2
double CMf = (ep - em) / (ep + 2 * em);
double acceptedtrials = 0., totaltrials = 0., acceptancerate;
int jmax = 200000000;
double rad, A, B, MetropolisMC, kBoltzmann = 1.38064852e-23, Temp = 293.;
double Diameter_post = 20., post_distance = 30.; //diameter of each post, tilted distance between posts  
double dist_posts_offset = 0, X_cyl1_offset = 18, Y_cyl1_offset = 18;
double number_of_slides_in_comsol_output = 30;

double post_H = 4, pH = 9.43, R = Di * (1.0e-6) / 2., eps = 78.5, e0 = 1.60217657e-19, Na = 6.022e+23, zetasphere = -43.88, zetapost = -23.88;
double A1w2, Asio2 = 6.6e-20, Aw = 3.7e-20;

int const Particle_num = 2500;
double paticles_cut_off_distance = 15.;

string input = "30zSlices-albanie-2post-1MHz-z=30-v=5-normal-results-mesh0.6-0.8.txt"; //input file name
string output = "out sio2 2500 particles di 20 dist 30 v5 1MHz.txt";//output file name
bool plot_dep_enrgy = false;//which type output we want
bool Xvalues = false, Yvalues = false, Zvalues = false, DEPvalues = true;
double Zd_cross = 6;
//const int numvalue = 565820;//as we know our data has this much row in it 565820


struct triple {

	vector<double> xvalue;
	vector<double> yvalue;
	vector<double> zvalue;
	vector<double> Enorm;
	//int totalnum; //this is only added to save the total number of data for being used out of function
};
triple getdata() {
	triple fromprevious;
	vector<double> values;
	double in = 0.0;

	ifstream infile;
	infile.open(input);//remember to save a text file in exact location that "out" file exist  
	while (infile >> in) {
		values.push_back(in);
	}
	int jj;// so imp you should update nn each time (which is the number points in each rectangle) NN is the number of particles
	for (jj = 0; jj < int(values.size()); jj++) {
		//cout << values[j] << endl;
		if (jj % 4 == 0) {
			fromprevious.xvalue.push_back(values[jj]);
			//fromprevious.xvalue[jj / 4] = values[jj];
		}
		if (jj % 4 == 1) {
			fromprevious.yvalue.push_back(values[jj]);
			//fromprevious.yvalue[jj / 4] = values[jj];
		}
		if (jj % 4 == 2) {
			fromprevious.zvalue.push_back(values[jj]);
			//fromprevious.zvalue[jj / 4] = values[jj];// it rotates periodiacally from x-y-z
		}
		if (jj % 4 == 3) {
			fromprevious.Enorm.push_back(values[jj]);
			//fromprevious.Enorm[jj / 4] = values[jj];// E:V/m value is savd here
		}
	}
	//fromprevious.totalnum = jj;

	return fromprevious;
}

//i made this function because rand is useful for int variables, RAND_MAX is known in library of <stdlib.h>
double fRand(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}
// exactly here we can check the risk of overlap for rectangles with PBC situation, then in overlapchecker we need to do some modifications to reach our goal in PBC
double calcDistance(double x1, double x2, double y1, double y2, double z1, double z2, double lx, double ly)
{
	//nwe added: previous problem was i didn't add abs
	// no one has overlap problem
	if (abs(x1 - x2) <= (lx / 2) && abs(y1 - y2) <= (ly / 2)) {
		return sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2) + pow((z1 - z2), 2));
	}//x has overlap
	else if (abs(x1 - x2) >= (lx / 2) && abs(y1 - y2) <= (ly / 2)) {
		return sqrt(pow((lx - abs(x1 - x2)), 2) + pow((y1 - y2), 2) + pow((z1 - z2), 2));
	}//x,y has overlap
	else if (abs(x1 - x2) >= (lx / 2) && abs(y1 - y2) >= (ly / 2)) {
		return sqrt(pow((lx - abs(x1 - x2)), 2) + pow((ly - abs(y1 - y2)), 2) + pow((z1 - z2), 2));
	}//y,z has overlap
	else if (abs(x1 - x2) <= (lx / 2) && abs(y1 - y2) >= (ly / 2)) {
		return sqrt(pow((x1 - x2), 2) + pow((ly - abs(y1 - y2)), 2) + pow((z1 - z2), 2));
	}
}
double Upfde(double Es)
{
	double landa;// , EE0;
	//EE0 = pow(8., -0.5)*Vpp / dg;

	landa = pi * em*pow((Di*1e-6 / 2.), 3.)*pow(CMf, 2.) / kT;
	return -2. * kT*landa*pow(Es, 2.) / CMf;
	//return 2. * pi * em*CMf *pow((Di*1e-6 / 2.), 3.)*pow(Es, 2.);

}
struct xy
{
	double x;
	double y;

};
class Energy_extraction {
public:// for error

	//Energy_extraction(float tempS) {
	//	sd = tempS;
	//	//Default=circleSize;
	//	generator = new Random();
	//}

	xy coordinate_transform(double Wcell_input, double Wcell_large, double x, double y)
	{

		double Di_p = Diameter_post, d = post_distance; //diameter of each post, tilted distance between posts  
		//vrtical or horizental distance between two posts
		double dist_posts = dist_posts_offset + (Di_p + d) / pow(2, 0.5);
		//X_cyl1 is th x coordinate of the bottom and left post in the raw data, same for Y_cyl1
		double X_cyl1 = X_cyl1_offset + Wcell_input / 2. - dist_posts, Y_cyl1 = Y_cyl1_offset + Wcell_input / 2. - dist_posts;

		//y is y in our large map, y_raw is what we have from data text file 
		//cell 1, 2,3,4
		if (y >= 0 && y < Wcell_large / 4.) {

			//cell 1
			if (x >= 0 && x < Wcell_large / 4.) {
				x = x + X_cyl1;
				y = y + Y_cyl1;
			}//cell 2
			else if (x >= Wcell_large / 4. && x < Wcell_large / 2.) {
				x = -x; x = x + Wcell_large / 2. + X_cyl1;
				y = y + Y_cyl1;

			}//cell 3
			else if (x >= Wcell_large / 2. && x < 3.*Wcell_large / 4.) {
				x = x - Wcell_large / 2. + X_cyl1;
				y = y + Y_cyl1;

			}//cell 4
			else if (x >= 3.*Wcell_large / 4. && x <= Wcell_large) {
				x = -x; x = x + Wcell_large + X_cyl1;
				y = y + Y_cyl1;

			}
		}
		//cells 5, 6,7,8
		else if (y >= Wcell_large / 4. && y < Wcell_large / 2.) {

			//cell 5
			if (x >= 0 && x < Wcell_large / 4.) {
				x = -x; x = x + Wcell_large / 4. + X_cyl1;
				y = y - Wcell_large / 4. + Y_cyl1;
			}//cell 6
			else if (x >= Wcell_large / 4. && x < Wcell_large / 2.) {
				x = x - Wcell_large / 4. + X_cyl1;
				y = y - Wcell_large / 4. + Y_cyl1;
			}//cell 7
			else if (x >= Wcell_large / 2. && x < 3.*Wcell_large / 4.) {
				x = -x; x = x + 3 * Wcell_large / 4. + X_cyl1;
				y = y - Wcell_large / 4. + Y_cyl1;
			}//cell 8
			else if (x >= 3.*Wcell_large / 4. && x <= Wcell_large) {
				x = x - 3 * Wcell_large / 4. + X_cyl1;
				y = y - Wcell_large / 4. + Y_cyl1;

			}
		}
		//cell 9, 10,11,12
		else if (y >= Wcell_large / 2. && y < 3 * Wcell_large / 4.) {
			//cell 9
			if (x >= 0 && x < Wcell_large / 4.) {
				x = x + X_cyl1;
				y = y - Wcell_large / 2. + Y_cyl1;
			}//cell 10
			else if (x >= Wcell_large / 4. && x < Wcell_large / 2.) {
				x = -x; x = x + Wcell_large / 2. + X_cyl1;
				y = y - Wcell_large / 2. + Y_cyl1;

			}//cell 11
			else if (x >= Wcell_large / 2. && x < 3.*Wcell_large / 4.) {
				x = x - Wcell_large / 2. + X_cyl1;
				y = y - Wcell_large / 2. + Y_cyl1;

			}//cell 12
			else if (x >= 3.*Wcell_large / 4. && x <= Wcell_large) {
				x = -x; x = x + Wcell_large + X_cyl1;
				y = y - Wcell_large / 2. + Y_cyl1;

			}
		}
		//cells 13, 14,15,16
		else if (y >= 3 * Wcell_large / 4. && y < Wcell_large) {

			//cell 13
			if (x >= 0 && x < Wcell_large / 4.) {
				x = -x; x = x + Wcell_large / 4. + X_cyl1;
				y = y - 3 * Wcell_large / 4. + Y_cyl1;
			}//cell 14
			else if (x >= Wcell_large / 4. && x < Wcell_large / 2.) {
				x = x - Wcell_large / 4. + X_cyl1;
				y = y - 3 * Wcell_large / 4. + Y_cyl1;
			}//cell 15
			else if (x >= Wcell_large / 2. && x < 3.*Wcell_large / 4.) {
				x = -x; x = x + 3 * Wcell_large / 4. + X_cyl1;
				y = y - 3 * Wcell_large / 4. + Y_cyl1;
			}//cell 16
			else if (x >= 3.*Wcell_large / 4. && x <= Wcell_large) {
				x = x - 3 * Wcell_large / 4. + X_cyl1;
				y = y - 3 * Wcell_large / 4. + Y_cyl1;

			}
		}
		xy result = { x, y };
		return result;
	}


	double Energy_Value_Extaction(double DataNum, double DataSection, double th, double Xd, double Yd, double Zd, vector<double> X, vector<double> Y, vector<double> Z, vector<double> E, double data_lx, double data_ly, double data_lz, double d_lx)
	{
		//in here we relate the value of Xd and Yd to our data coordinate
		Xd = coordinate_transform(data_lx, d_lx, Xd, Yd).x;
		Yd = coordinate_transform(data_lx, d_lx, Xd, Yd).y;

		//th=threshold
		double xdd, Xdd, xdd_, Xdd_, ydd, Ydd, ydd_, Ydd_, zdd, Zdd_;//Xd=Xdesird, Yd=Ydesired, Zd=Zdesired
		int mm = 0; double Ed = 0.0;//Ed=Ederivd

		int i = 0, o = 0, imin, kmin, imax, kmax;
		xdd = Xd - th, Xdd = Xd, xdd_ = Xd, Xdd_ = Xd + th; //Ydd = Yd - th, Ydd_ = Yd + th; Zdd = Zd - 2., Zdd_ = Zd + 2.;
		ydd = Yd - th, Ydd = Yd, ydd_ = Yd, Ydd_ = Yd + th;
		zdd = -0.0001, Zdd_ = Zd + 0.0001;//some points have values only for z=0 or lz=8, so we need to cover all the range from 0 to lz

		if (xdd < 0.) { xdd = xdd + data_lx, Xdd = data_lx + 0.0001, xdd_ = -0.0001; } if (ydd < 0.) { ydd = ydd + data_ly, Ydd = data_ly + 0.0001, ydd_ = -0.0001; }// if (zdd < 0.) { zdd = -0.0001; };
		if (Xdd_ > data_lx) { xdd_ = -0.0001, Xdd = data_lx + 0.0001, Xdd_ = Xdd_ - data_lx; } if (Ydd_ > data_ly) { ydd_ = -0.0001, Ydd = data_ly + 0.0001, Ydd_ = Ydd_ - data_ly; }// if (Zdd_ > lz) { Zdd_ = lz + 0.0001; };

		double A1 = 100, A2 = 100, A3 = 100;//A1 and A2 are the colsest and the second closest distance from the point (Xd, Yd, Zd), n1 and n2 are the index of these points

		int n1 = 1, n2 = 2, n3 = 3;
		imax = DataNum, imin = 0;

		//here I have made the code compatable for the comsol slice data

		double number_of_slides = number_of_slides_in_comsol_output;
		double slide_space = round(DataNum / number_of_slides);
		//Zd = 30;
		double gh = round(Zd);
		if (gh == round(5)) {
			gh = round(6);
		}
		//int DataSection2 = int(slide_space / data_lx);//140062/70.71=1981  |3781=567224/150, 565806 are the total # of data and we know x coordinate is sorted, by doing this we can make the code faster

		//double mnnm = gh* slide_space+DataSection2;
		//double bbb = Z.at(mnnm);

		//double mnnm1 = (gh-1) * slide_space + DataSection2 * int(Xd - th - 3);
		//double bbb2 = Z.at(mnnm1);

		//double mnnm2 = gh * slide_space + DataSection2 * int(Xd + th + 3);
		//double bbb3 = Z.at(mnnm2), bbb4 = Z.at(slide_space), bbb5 = Z.at(2*slide_space), bbb6 = Z.at(3 * slide_space), bbb7 = Z.at(4 * slide_space), bbb8 = Z.at(5 * slide_space), bbb9 = Z.at(6 * slide_space);
		//double bbb10 = Z.at(7*slide_space), bbb11 = Z.at(8 * slide_space), bbb12 = Z.at(9 * slide_space), bbb13 = Z.at(10 * slide_space), bbb14 = Z.at(11 * slide_space), bbb15 = Z.at(12 * slide_space);
		//double bbb16 = Z.at(13*slide_space), bbb17 = Z.at(14 * slide_space), bbb18 = Z.at(15 * slide_space), bbb19 = Z.at(16 * slide_space), bbb20 = Z.at(17 * slide_space), bbb21 = Z.at(18* slide_space);
		//double bbb00 = X.at(4 * slide_space);

		if ((Zd - th) > 0 && (Zd + th) < (data_lz)) {
			imin = (gh - th + 1) * slide_space, imax = (gh + th - 2) * slide_space;//140062/70.71=1981  |3781=567224/150, 565806 are the total # of data and we know x coordinate is sorted, by doing this we can make the code faster
		}
		else if (Zd <= th) {
			imax = (gh + th - 2) * slide_space;
		}
		else if (Zd >= (data_lz - th)) {
			imin = (gh - th + 1) * slide_space;
		}
		//double bbb00 = Z.at(imin);
		//double bbb01 = Z.at(imax-1);

		for (i = imin; i < imax; i++) {

			//if ((round(Z.at(i)) == gh) ) {

			if ((X.at(i) >= xdd && X.at(i) <= Xdd) || (X.at(i) >= xdd_ && X.at(i) <= Xdd_))


			{

				if ((Y.at(i) >= ydd && Y.at(i) <= Ydd) || (Y.at(i) >= ydd_ && Y.at(i) <= Ydd_))
				{


					//if (Z.at(i) > z0 && Z.at(i) <= (data_lz + 0.001))
					//{
					mm++;

					if (mm == 1) { A1 = calcDistance(Xd, X.at(i), Yd, Y.at(i), Zd, Z.at(i), data_lx, data_ly), n1 = i; }
					if (mm == 2) {
						A2 = calcDistance(Xd, X.at(i), Yd, Y.at(i), Zd, Z.at(i), data_lx, data_ly), n2 = i;;
						if (A2 < A1) { A3 = A1, A1 = A2, A2 = A3, n3 = n1, n1 = n2, n2 = n3; }//now we know for sur that A1<A2
					}
					if (mm > 2) { //here we want to updat A1 and A2 to mak sure we'v got the closest points
						A3 = calcDistance(Xd, X.at(i), Yd, Y.at(i), Zd, Z.at(i), data_lx, data_ly), n3 = i;

						if (A3 < A1) { A2 = A1, A1 = A3, n2 = n1, n1 = n3; }
						else if ((A3 < A2) && (A3 >= A1)) { A2 = A3, n2 = n3; }
					}

					//}
				}
			}
			//}

		}

		return Ed = (E.at(n1) + E.at(n2)) / 2.;


	}





};



class mutual_interactions {
public:
	double ys(double zeta)
	{
		return (zeta*e0) / (kT);
	}
	double kappa(double pH)
	{
		double c = pow(10., -pH)* 1.e+3;
		return sqrt((2.* pow(e0, 2.)*Na*c) / (eps* eps0* kT));
	}
	double yy(double r, double zetaau)
	{
		return 4.*exp(0.5 *kappa(pH)* (r - 2.*R)) * atanh(exp(-1.*kappa(pH)*(r - 2.*R) / 2.)* tanh(ys(zetaau) / 4.));

	}
	//electroestatic between rod-seg rod-seg
	//double ESt_SR_SR(double r, double zeta1, double zeta2)
	double ESt_SR_SR(double r)
	{
		r = r * 1.0e-6;
		double zeta1 = zetasphere, zeta2 = zetasphere;
		return eps * eps0*pow((kT / e0), 2.)* yy(r, zeta1)*yy(r, zeta2)*(pow(R, 2.) / r)*log(1 + exp(-kappa(pH) * (r - 2.*R)));//modified eq for wo sphrs with differnt zeta potntial (particle spheres and post spheres) 

	}
	double Hr(double x, double y)
	{

		return y / (pow(x, 2.) + x * y + x) + y / (pow(x, 2.) + x * y + x + y) + 2.*log((pow(x, 2.) + x * y + x) / (pow(x, 2.) + x * y + x + y));
	}
	//double EvdwSS(double r, double Asphere1, double Asphere2)
	double EvdwSS(double r)
	{//it includes all core shell, shell shell, core core interactions
		r = r * 1.0e-6;
		double Asphere1, Asphere2;
		Asphere1 = Asio2;
		Asphere2 = Asio2;
		A1w2 = (pow(Asphere1, 1. / 2) - pow(Aw, 1. / 2)) * (pow(Asphere2, 1. / 2) - pow(Aw, 1. / 2));
		//no problem if R and ds are not in nm scale because x and y in Hr function are ratios
		//encourted a strange problem don't write 1/12 will make zero in results, use 1./12 !!!!
		return (-1. / 12)* (A1w2*Hr((r - 2.*R) / (2.*R), 1.));
	}
	double EDL_S_FP_bevan(double h)//Electric double layer between Sphere and Flat Plate Equation
	{
		double zeta1 = zetasphere, zeta2 = zetapost;//for now we consider the post zeta as -23.88
		double sai1 = tanh(ys(zeta1) / 4.);
		double sai2 = tanh(ys(zeta2) / 4.);
		return 64 * pi *R* eps*eps0*sai1*sai2*pow((kT / e0), 2.)*exp(-kappa(pH) * ((h - post_H)* (1.0e-6) - R));//eq from SEI paper: Bhattacharjee, Elimelech
		//return r;
	}
	double EDL_S_FP(double h)//Electric double layer between Sphere and Flat Plate Equation
	{
		double zeta1 = zetasphere, zeta2 = zetapost;//for now we consider the post zeta as -23.88
		double sai1 = tanh(ys(zeta1) / 4.);
		double sai2 = tanh(ys(zeta2) / 4.);
		return 32 * eps*eps0*kappa(pH)*sai1*sai2*pow((kT / e0), 2.)*exp(-kappa(pH) * (h));//eq from SEI paper: Bhattacharjee, Elimelech
		//return r;
	}
	//a function to calculate integral numarically
	double SEIntegralSumEq(double lowBound, double upBound, int n, double z)//Surface Element Integration Equation also the integral is calculated numarically 
	{
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
			double funcValue = EDL_S_FP(h1) - EDL_S_FP(h2);
			double rectangleArea = funcValue * r*dr;
			cumSum += rectangleArea;

		}
		double yyy = 2 * pi*cumSum;
		return 2 * pi*cumSum;
	}

};


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
	ofstream out(output);
	//streambuf *coutbuf = std::cout.rdbuf();
	cout.rdbuf(out.rdbuf());


	//making the initial position of particles
	srand(time(NULL));


	triple xyz;
	xyz = getdata();
	//int totnum, xvalnum;
	//totnum = xyz.totalnum;
	//xvalnum = int(totnum/4)+10;
	//N = xvalnum;
	Energy_extraction energy_extraction;
	mutual_interactions mutual_interactions;
	//for (it = X.begin(); it != X.end(); ++it)	{	
		//X.emplace(it + 5, xyz.xvalue[it]);
		//X[i] = xyz.xvalue[i];
		//Y[i] = xyz.yvalue[i];
		//Z[i] = xyz.zvalue[i];
		//E[i] = xyz.Enorm[i];
	//}
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




	int const Pnum = Particle_num;
	int k = 0, ke = 0;
	double Px[Pnum][2], Py[Pnum][2], Pz[Pnum][2], PEnergy[Pnum][2], A, B, A1, A2, A3, A4, A5, A6;//Di is in 1um
	int OverlapChance, ff = 0;
	double G;
	if (plot_dep_enrgy == false) {
		for (i = 0; i < Pnum; i++) {
			Px[i][0] = (double)fRand(0., d_lx);
			Py[i][0] = (double)fRand(0., d_ly);

			//debugged: in direction of Z box will start from 0 upto lz not from -lz/2 upto lz/2
			//new debugged: based on exp(-kappa * (z - (Di / 2.))) z cannot be lower than radius of particle 'cause it will pass the beneath wall which is not real
			Pz[i][0] = (double)fRand(z0 + 0.0001, d_lz);



			for (k = (i - 1); k >= 0; k--) {

				int OverlapChance;
				A = calcDistance(Px[i][0], Px[k][0], Py[i][0], Py[k][0], Pz[i][0], Pz[k][0], d_lx, d_ly);
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

		for (k = 0; k < Pnum; k++) {
			G = (double)zero;
			//first we need to extract enrgy value of each point from our text data
			Xd = Px[k][0], Yd = Py[k][0], Zd = Pz[k][0];
			Ed = energy_extraction.Energy_Value_Extaction(DataNum, DataSection, th, Xd, Yd, Zd, X, Y, Z, E, data_lx, data_ly, data_lz, d_lx);
			//now we add Electrostatic among particles
			for (ke = 0; ke < Pnum; ke++) {
				if (ke == k) {
					ke++;
				}
				A = calcDistance(Px[k][0], Px[ke][0], Py[k][0], Py[ke][0], Pz[k][0], Pz[ke][0], d_lx, d_ly);
				B = paticles_cut_off_distance;
				// new added: we need to add A > 340 because if there's a overlab the energy of system will increase to infinity						
				if (A <= B) {
					//G += EvdwSS(A) + ESt_SR_SR(A) + ESt_SR_ER(A) + ESt_ER_ER(A);
					//G[i][1] += EvdwSS(A) + ESt_SR_SR(A);								
					G += mutual_interactions.EvdwSS(A) + mutual_interactions.ESt_SR_SR(A);
				}
			}
			A1 = calcDistance(Px[k][0], (d_lx / 2.) - 0, Py[k][0], (d_ly / 2.) - 0, Pz[k][0], Pz[k][0], d_lx, d_ly);//we consider all z values equal as we only care about 2D distance of particles from the center of the posts  
			A2 = calcDistance(Px[k][0], (d_lx / 2.) - bb, Py[k][0], (d_ly / 2.) - bb, Pz[k][0], Pz[k][0], d_lx, d_ly);
			A3 = calcDistance(Px[k][0], (d_lx / 2.) - bb, Py[k][0], (d_ly / 2.) + bb, Pz[k][0], Pz[k][0], d_lx, d_ly);
			A4 = calcDistance(Px[k][0], (d_lx / 2.) + bb, Py[k][0], (d_ly / 2.) - bb, Pz[k][0], Pz[k][0], d_lx, d_ly);
			A5 = calcDistance(Px[k][0], (d_lx / 2.) + bb, Py[k][0], (d_ly / 2.) + bb, Pz[k][0], Pz[k][0], d_lx, d_ly);
			A6 = calcDistance(Px[k][0], 0, Py[k][0], 0, Pz[k][0], Pz[k][0], d_lx, d_ly);

			if ((A1 <= Diameter / 2.) || (A2 <= Diameter / 2.) || (A3 <= Diameter / 2.) || (A4 <= Diameter / 2.) || (A5 <= Diameter / 2.) || (A6 <= Diameter / 2.)) {
				G += mutual_interactions.EDL_S_FP_bevan(Zd);
			}

			PEnergy[k][0] = Upfde(Ed) + G;
		}

		int j = 0;
		// NU is a counter to track usage of norm[] array and if it end, so we renew it
		for (j = 0; j < jmax; j++) {
			totaltrials++;
			int min = 0, max = Pnum - 1;

			int randi = rand() % (max - min + 1) + min;

			double radnum, randdeltax, randdeltay, randdeltaz;
			radnum = (double)fRand(0., 1.);
			if (0. <= radnum && radnum < (1. / 3)) {

				randdeltax = (double)fRand(-d_lx / 2.5, d_lx / 2.5);
				Px[randi][1] = Px[randi][0] + randdeltax;
				Py[randi][1] = Py[randi][0];
				Pz[randi][1] = Pz[randi][0];

			}
			else if ((1. / 3) <= radnum && radnum < (2. / 3)) {

				randdeltay = (double)fRand(-d_ly / 2.5, d_ly / 2.5);
				Px[randi][1] = Px[randi][0];
				Py[randi][1] = Py[randi][0] + randdeltay;
				Pz[randi][1] = Pz[randi][0];

			}
			else {
				randdeltaz = (double)fRand(-(d_lz - z0) / 5., (d_lz - z0) / 5.);
				Px[randi][1] = Px[randi][0];
				Py[randi][1] = Py[randi][0];
				Pz[randi][1] = Pz[randi][0] + randdeltaz;
			}

			//pbc
			//pbc for x direction
			if (Px[randi][1] > d_lx) {

				Px[randi][1] = Px[randi][1] - d_lx;
			}
			else if (Px[randi][1] < 0.) {

				Px[randi][1] = Px[randi][1] + d_lx;
			}
			//pbc for y direction
			if (Py[randi][1] > d_ly) {

				Py[randi][1] = Py[randi][1] - d_ly;
			}
			else if (Py[randi][1] < 0.) {

				Py[randi][1] = Py[randi][1] + d_ly;
			}

			//pbc for z direction
			if (Pz[randi][1] > d_lz) {

				//Pz[randi][1] = lz;
				Pz[randi][1] = Pz[randi][1] - (d_lz - z0);
			}
			else if (Pz[randi][1] < z0) {
				//this should not happen as we have always reulsion with the below surface 
				//Pz[randi][1] = z0+0.001;
				Pz[randi][1] = Pz[randi][1] + (d_lz - z0);

			}
			o = 0;
			for (i = 0; i < Pnum; i++) {
				ff = 0;
				if (i != randi) {
					A = calcDistance(Px[randi][1], Px[i][0], Py[randi][1], Py[i][0], Pz[randi][1], Pz[i][0], d_lx, d_ly);
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
						Xd = Px[randi][1], Yd = Py[randi][1], Zd = Pz[randi][1];
						//now a loop to find PEnergy value by searching through the Edelta2 values 
						Ed = energy_extraction.Energy_Value_Extaction(DataNum, DataSection, th, Xd, Yd, Zd, X, Y, Z, E, data_lx, data_ly, data_lz, d_lx);
						for (ke = 0; ke < Pnum; ke++) {
							if (ke != randi) {
								A = calcDistance(Px[randi][1], Px[ke][1], Py[randi][1], Py[ke][1], Pz[randi][1], Pz[ke][1], d_lx, d_ly);
								B = paticles_cut_off_distance;
								// new added: we need to add A > 340 because if there's a overlab the energy of system will increase to infinity						
								if (A <= B) {
									//G += EvdwSS(A) + ESt_SR_SR(A) + ESt_SR_ER(A) + ESt_ER_ER(A);
									//G[i][1] += EvdwSS(A) + ESt_SR_SR(A);								
									G += mutual_interactions.ESt_SR_SR(A) + mutual_interactions.EvdwSS(A);

								}
							}
						}
						A1 = calcDistance(Px[randi][1], (d_lx / 2.) - 0, Py[randi][1], (d_ly / 2.) - 0, Pz[randi][1], Pz[randi][1], d_lx, d_ly);
						A2 = calcDistance(Px[randi][1], (d_lx / 2.) - bb, Py[randi][1], (d_ly / 2.) - bb, Pz[randi][1], Pz[randi][1], d_lx, d_ly);
						A3 = calcDistance(Px[randi][1], (d_lx / 2.) - bb, Py[randi][1], (d_ly / 2.) + bb, Pz[randi][1], Pz[randi][1], d_lx, d_ly);
						A4 = calcDistance(Px[randi][1], (d_lx / 2.) + bb, Py[randi][1], (d_ly / 2.) - bb, Pz[randi][1], Pz[randi][1], d_lx, d_ly);
						A5 = calcDistance(Px[randi][1], (d_lx / 2.) + bb, Py[randi][1], (d_ly / 2.) + bb, Pz[randi][1], Pz[randi][1], d_lx, d_ly);
						A6 = calcDistance(Px[randi][1], 0, Py[randi][1], 0, Pz[randi][1], Pz[randi][1], d_lx, d_ly);

						if ((A1 <= Diameter / 2.) || (A2 <= Diameter / 2.) || (A3 <= Diameter / 2.) || (A4 <= Diameter / 2.) || (A5 <= Diameter / 2.) || (A6 <= Diameter / 2.)) {
							//double Gu= mutual_interactions.SEIntegralSumEq(0, R, 500, Zd);
							//double Guu = mutual_interactions.EDL_S_FP_bevan(Zd);
							G += mutual_interactions.EDL_S_FP_bevan(Zd);
						}

						PEnergy[randi][1] = Upfde(Ed) + G;

						if (PEnergy[randi][0] >= PEnergy[randi][1]) {
							// false == 0 and true = !false
							MetropolisMC = true;
						}

						else {
							double P = exp(-(PEnergy[randi][1] - PEnergy[randi][0]) / (kBoltzmann*Temp));
							//double P = exp(-(TotalEnergy[1] - TotalEnergy[0]) / (kBoltzmann*Temp));
							rad = (double)fRand(0., 1.);
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

							Px[randi][0] = Px[randi][1];
							Py[randi][0] = Py[randi][1];
							Pz[randi][0] = Pz[randi][1];
							PEnergy[randi][0] = PEnergy[randi][1];



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
									cout << "{RGBColor[224, 255, 255], Opacity[0.8], Sphere[{" << Px[i][0] << ", " << Py[i][0] << ", " << Pz[i][0] << "}, " << (Di / 2.) << "]},";
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
			}

		}
	}
	if (plot_dep_enrgy == true) {
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
					cout << Upfde(Ed) << "\n"; ;
				}
			}
		}
	}
	return 0;
}
