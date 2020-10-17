#include "Energy_extraction.h"
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


MC_funcs MC;//calling the MC_funcs class inside our Energy_extraction class as a private element

Energy_extraction::xy Energy_extraction::coordinate_transform(double Wcell_input, double Wcell_large, double x, double y) {

	double Di_p = Diameter_post, d = post_distance; //diameter of each post, tilted distance between posts  
	//vrtical or horizental distance between two posts
	double dist_posts = (Di_p + d) / pow(2, 0.5);
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
	Energy_extraction::xy result = { x, y };
	return result;
}
	


double Energy_extraction::Energy_Value_Extaction(double DataNum, double DataSection, double th, double Xd, double Yd, double Zd, vector<double> X, vector<double> Y, vector<double> Z, vector<double> E, double data_lx, double data_ly, double data_lz, double d_lx) {
	//in here we relate the value of Xd and Yd to our data coordinate
	Xd = Energy_extraction::coordinate_transform(data_lx, d_lx, Xd, Yd).x;
	Yd = Energy_extraction::coordinate_transform(data_lx, d_lx, Xd, Yd).y;

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

				if (mm == 1) { A1 = MC.calcDistance(Xd, X.at(i), Yd, Y.at(i), Zd, Z.at(i), data_lx, data_ly), n1 = i; }
				if (mm == 2) {
					A2 = MC.calcDistance(Xd, X.at(i), Yd, Y.at(i), Zd, Z.at(i), data_lx, data_ly), n2 = i;;
					if (A2 < A1) { A3 = A1, A1 = A2, A2 = A3, n3 = n1, n1 = n2, n2 = n3; }//now we know for sur that A1<A2
				}
				if (mm > 2) { //here we want to updat A1 and A2 to mak sure we'v got the closest points
					A3 = MC.calcDistance(Xd, X.at(i), Yd, Y.at(i), Zd, Z.at(i), data_lx, data_ly), n3 = i;

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


