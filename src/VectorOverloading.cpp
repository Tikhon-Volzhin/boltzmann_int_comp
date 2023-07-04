#include <iostream>
#include <vector>
#include <random>
#include <string> 
#include <fstream>
#include <cmath>

#include "headers/VectorOverloading.h"
#include "headers/AdditionalFunctions.h"
#define Vec_3d std::vector<vector<vector<double>>>

using namespace std;


vector<double> operator + (const vector<double>& v1, const vector<double>& v2)
{
	const size_t len_max = min<size_t>(v1.size(), v2.size());
	vector<double> result(len_max);

	for (size_t i = 0; i < result.size() && i < v2.size(); i++)
	{
		result[i] = v1[i] + v2[i];
	}

	return result;
}
vector<int> operator + (const vector<int>& v1, const vector<int>& v2)
{
	const size_t len_max = min<size_t>(v1.size(), v2.size());
	vector<int> result(len_max);

	for (size_t i = 0; i < result.size() && i < v2.size(); i++)
	{
		result[i] = v1[i] + v2[i];
	}

	return result;
}


vector<double> operator - (const vector<double>& v1, const vector<double>& v2)
{
	const size_t len_max = min<size_t>(v1.size(), v2.size());
	vector<double> result(len_max);

	for (size_t i = 0; i < result.size() && i < v2.size(); i++)
	{
		result[i] = v1[i] - v2[i];
	}

	return result;
}
vector<int> operator - (const vector<int>& v1, const vector<int>& v2)
{
	const size_t len_max = min<size_t>(v1.size(), v2.size());
	vector<int> result(len_max);

	for (size_t i = 0; i < result.size() && i < v2.size(); i++)
	{
		result[i] = v1[i] - v2[i];
	}

	return result;
}


double operator ^ (const vector<double>& v1, const vector<double>& v2)
{
	double result = 0;

	for (size_t i = 0; i < v1.size() && i < v2.size(); i++)
	{
		result += v1[i] * v2[i];
	}

	return result;
}
int operator ^ (const vector<int>& v1, const vector<int>& v2)
{
	int result = 0;

	for (size_t i = 0; i < v1.size() && i < v2.size(); i++)
	{
		result += v1[i] * v2[i];
	}

	return result;
}


vector<double> operator * (const vector<double>& v1, const vector<double>& v2)
{
	const size_t len_max = min<size_t>(v1.size(), v2.size());
	vector<double> result(len_max);

	for (size_t i = 0; i < result.size() && i < v2.size(); i++)
	{
		result[i] = v1[i]*v2[i];
	}

	return result;
}
vector<int> operator * (const vector<int>& v1, const vector<int>& v2)
{
	const size_t len_max = min<size_t>(v1.size(), v2.size());
	vector<int> result(len_max);

	for (size_t i = 0; i < result.size() && i < v2.size(); i++)
	{
		result[i] = v1[i] * v2[i];
	}

	return result;
}
vector<double> operator * (const double real_number, const vector<double>& vec)
{
	vector<double> result = vec;

	for (size_t i = 0; i < vec.size(); i++)
	{
		result[i] *= real_number;
	}

	return result;
}

vector<double> operator * (const double real_number, const vector<int>& vec)
{
	vector<double> result(vec.size());

	for (size_t i = 0; i < vec.size(); i++)
	{
		result[i] = real_number* (double)vec[i];
	}

	return result;
}
vector<int> operator * (const int real_number, const vector<int>& vec)
{
	vector<int> result(vec.size());

	for (size_t i = 0; i < vec.size(); i++)
	{
		result[i] = real_number * vec[i];
	}

	return result;
}



vector<double> operator / (const vector<double>& vec, const double real_number)
{
	vector<double> result(vec.size());

	for (size_t i = 0; i < vec.size(); i++)
	{
		result[i] = vec[i] / real_number;
	}

	return result;
}
vector<double> operator / (const vector<int>& vec, const double real_number)
{
	vector<double> result(vec.size());

	for (size_t i = 0; i < vec.size(); i++)
	{
		result[i] = (double)vec[i]/real_number;
	}

	return result;
}


bool operator ==(const vector<double>& v1, const vector<double>& v2)
{
	if (v1.size() != v2.size())
		return false;

	for (size_t i = 0; i < v1.size(); i++)
	{
		if (!are_same_db_weak(v1[i], v2[i]))
			return false;
	}
	return true;
}
bool operator ==(const vector<int>& v1, const vector<int>& v2)
{
	if (v1.size() != v2.size()) 
		return false;
	
	for (size_t i = 0; i < v1.size(); i++)
	{
		if (v1[i] != v2[i]) 
			return false;
	}
	return true;
}



Vec_3d operator + (const Vec_3d& v1, const Vec_3d& v2) {
	int N_x = v1.size(), N_y = v1[0].size(), N_z = v1[0][0].size();
	Vec_3d result(N_x, vector<vector<double>>(N_y, vector<double>(N_z)));
	if (N_x == v2.size() && N_y == v2[0].size() && N_z == v2[0][0].size()) {
		for (int i = 0; i < N_x; i++) {
			for (int j = 0; j < N_y; j++) {
				for (int k = 0; k < N_z; k++) {
					result[i][j][k] = (v1[i][j][k] + v2[i][j][k]);
				}
			}
		}
	}
	return result;
}
Vec_3d operator - (const Vec_3d& v1, const Vec_3d& v2) {
	int N_x = v1.size(), N_y = v1[0].size(), N_z = v1[0][0].size();
	Vec_3d result(N_x, vector<vector<double>>(N_y, vector<double>(N_z)));
	if (N_x == v2.size() && N_y == v2[0].size() && N_z == v2[0][0].size()) {
		for (int i = 0; i < N_x; i++) {
			for (int j = 0; j < N_y; j++) {
				for (int k = 0; k < N_z; k++) {
					result[i][j][k] = (v1[i][j][k] - v2[i][j][k]);
				}
			}
		}
	}
	return result;
}
Vec_3d operator / (const Vec_3d& v1, const Vec_3d& v2) {
	int N_x = v1.size(), N_y = v1[0].size(), N_z = v1[0][0].size();
	Vec_3d result(N_x, vector<vector<double>>(N_y, vector<double>(N_z)));
	if (N_x == v2.size() && N_y == v2[0].size() && N_z == v2[0][0].size()) {
		for (int i = 0; i < N_x; i++) {
			for (int j = 0; j < N_y; j++) {
				for (int k = 0; k < N_z; k++) {
					if (v2[i][j][k] == 0)
						result[i][j][k] = 0;
					else
						result[i][j][k] = (v1[i][j][k]/v2[i][j][k]);
				}
			}
		}
	}
	return result;
}

Vec_3d operator * (double scalar, const Vec_3d& v1) {
	int N_x = v1.size(), N_y = v1[0].size(), N_z = v1[0][0].size();
	Vec_3d result(N_x, vector<vector<double>>(N_y, vector<double>(N_z)));
	for (int i = 0; i < N_x; i++) {
		for (int j = 0; j < N_y; j++) {
			for (int k = 0; k < N_z; k++) {
				result[i][j][k] = v1[i][j][k] * scalar;
			}
		}
	}
	return result;
}
Vec_3d operator * (int scalar, const Vec_3d& v1) {
	int N_x = v1.size(), N_y = v1[0].size(), N_z = v1[0][0].size();
	Vec_3d result(N_x, vector<vector<double>>(N_y, vector<double>(N_z)));
	for (int i = 0; i < N_x; i++) {
		for (int j = 0; j < N_y; j++) {
			for (int k = 0; k < N_z; k++) {
				result[i][j][k] = v1[i][j][k] * (double)scalar;
			}
		}
	}
	return result;
}
Vec_3d operator / (const Vec_3d& v1, double scalar) {
	int N_x = v1.size(), N_y = v1[0].size(), N_z = v1[0][0].size();
	Vec_3d result(N_x, vector<vector<double>>(N_y, vector<double>(N_z)));
		for (int i = 0; i < N_x; i++) {
			for (int j = 0; j < N_y; j++) {
				for (int k = 0; k < N_z; k++) {
					result[i][j][k] = v1[i][j][k]/scalar;
				}
			}
		}
	return result;
}
Vec_3d operator / (const Vec_3d& v1, int scalar) {
	int N_x = v1.size(), N_y = v1[0].size(), N_z = v1[0][0].size();
	Vec_3d result(N_x, vector<vector<double>>(N_y, vector<double>(N_z)));
	for (int i = 0; i < N_x; i++) {
		for (int j = 0; j < N_y; j++) {
			for (int k = 0; k < N_z; k++) {
				result[i][j][k] = v1[i][j][k]/scalar;
			}
		}
	}
	return result;
}


double operator ^ (const Vec_3d& v1, const Vec_3d& v2) {
	int N_x = v1.size(), N_y = v1[0].size(), N_z = v1[0][0].size();
	double result = 0;
	if (N_x == v2.size() && N_y == v2[0].size() && N_z == v2[0][0].size()) {
		for (int i = 0; i < N_x; i++) {
			for (int j = 0; j < N_y; j++) {
				for (int k = 0; k < N_z; k++) {
					result += v1[i][j][k]* v2[i][j][k];
				}
			}
		}
	}
	return result;
}

