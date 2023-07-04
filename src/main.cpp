#include <iostream>
#include <vector>
#include <mutex>
#include <random>
#include <cmath>
#include <iomanip>

#include "headers/InitConditions.h"
#include "headers/VectorOverloading.h"
#include "headers/Net.h"
#include "headers/AdditionalFunctions.h"

double pi = 3.14159265358979323846;
double f_1_integral_factor = 1.0 / 4.0 / pi/ sqrt(2.0 * pi);
double f_2_integral_factor = 1.0 / 2.0 / pi/ sqrt(2.0 * pi);

double f_1(vector<double> & vec, double T, double dz_cut) {
	double squared_dz = vec^vec;
	
	if (squared_dz > dz_cut*dz_cut) return 0.0;
	return (exp(-0.5 * squared_dz) +  exp(-0.5 * squared_dz / T)/T / sqrt(T))* f_1_integral_factor;
}
double f_2(vector<double>& vec, double T, double dz_cut) {
	double squared_dz = vec ^ vec;
	
	if (squared_dz > dz_cut* dz_cut) return 0.0;
	if (vec[0] > 0) {
		return  exp(-0.5 * squared_dz)*f_2_integral_factor;
	}
	else {
		return  exp(-0.5 * squared_dz / T)/T/sqrt(T)* f_2_integral_factor;
	}
}

InitConditions F_1("F_1",f_1, 0.5, {2,2,2});
InitConditions F_2("F_2",f_2, 0.4, {1,2,2});


int main()
{
	Net net1(30, 1, F_2);
	net1.compute_num_I();
	return 0;
}
