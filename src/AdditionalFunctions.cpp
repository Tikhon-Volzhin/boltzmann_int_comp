#include <iostream>
#include <vector>
#include <random>
#include <string> 
#include <fstream>
#include <cmath>

#include "headers/AdditionalFunctions.h"


bool are_same_db_strict(double a, double b) {
	return abs(a - b) < 1e-30;
}
bool are_same_db_weak(double a, double b) {
	return abs(a - b) < 1e-12;
}
bool are_same_db_e5(double a, double b) {
	return abs(a - b) < 1e-12;
}
void write_3dvec_in_file(const Vec_3d& vec, string name, vector<int> cut) {
	ofstream file;
	file.open(name);
	int N_x = vec.size(), N_y = vec[0].size(), N_z = vec[0][0].size();
	file << N_x << " " << N_y << " " << N_z << endl;
	if (file.is_open()) {
		for (int i = 0; i < N_x / cut[0]; i++) {
			for (int j = 0; j < N_y / cut[1]; j++) {
				for (int k = 0; k < N_z / cut[2]; k++) {
					file << i << " " << j << " " << k << " " << vec[i][j][k] << endl;
				}
			}
		}
	}
	file.close();
}
Vec_3d read_3dvec_from_file(string name) {
	std::ifstream in(name);
	int N_x, N_y, N_z;
	in >> N_x >> N_y >> N_z;
	Vec_3d result(N_x, vector<vector<double>>(N_y, vector<double>(N_z)));
	int i, j, k;
	double f;
	if (in.is_open())
	{
		while (in >> i >> j >> k >> f)
		{
			result[i][j][k] = (double)f;
		}
	}
	in.close();
	return result;
}

std::ostream& operator<< (std::ostream& out, const vector<int>& vec)
{
	out << "[";
	for (int i = 0; i < vec.size() - 1; i++) {
		out << vec[i] << ", ";
	}
	out << vec[vec.size() - 1];
	out << "] ";
	return out;
}
std::ostream& operator<< (std::ostream& out, const vector<double>& vec)
{
	out << "[";
	for (int i = 0; i < vec.size() - 1; i++) {
		out << vec[i] << ", ";
	}
	out << vec[vec.size() - 1];
	out << "] ";
	return out;
}

vector<int> round_vec(const vector<double>& vec) {
	vector<int> result(vec.size());
	for (size_t i = 0; i < vec.size(); i++)
	{
		result[i] = (int)round(vec[i]);
	}
	return result;
}
vector<int> floor_vec(const vector<double>& vec) {
	vector<int> result(vec.size());
	for (size_t i = 0; i < vec.size(); i++)
	{
		result[i] = (int)floor(vec[i]);
	}
	return result;
}

double mean_3d_vec(const Vec_3d& vec) {
	int N_x = vec.size(), N_y = vec[0].size(), N_z = vec[0][0].size();
	double result = 0;
	int count = 0;
		for (int i = 0; i < N_x; i++) {
			for (int j = 0; j < N_y; j++) {
				for (int k = 0; k < N_z; k++) {
					if (vec[i][j][k] != 0) {
						result += vec[i][j][k];
						count++;
					}
				}
			}
		}
	return result/count;
}
double std_3d_vec(const Vec_3d& vec) {
	int N_x = vec.size(), N_y = vec[0].size(), N_z = vec[0][0].size();
	double result = 0;
	double mean = mean_3d_vec(vec);
	int count = 0;
	for (int i = 0; i < N_x; i++) {
		for (int j = 0; j < N_y; j++) {
			for (int k = 0; k < N_z; k++) {
				if (vec[i][j][k] != 0) {
					result += (vec[i][j][k] - mean) * (vec[i][j][k] - mean);
					count++;
				}
			}
		}
	}
	return sqrt(result / count);
}
double max_3d_vec(Vec_3d& vec) {
	int N_x = vec.size(), N_y = vec[0].size(), N_z = vec[0][0].size();
	double result = vec[0][0][0];
	for (int i = 0; i < N_x; i++) {
		for (int j = 0; j < N_y; j++) {
			for (int k = 0; k < N_z; k++) {
				if (vec[i][j][k] > result) {
					result = vec[i][j][k];
				}
			}
		}
	}
	return result;
}