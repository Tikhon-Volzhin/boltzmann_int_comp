#include <iostream>
#include <vector>

#include "headers/InitConditions.h"

InitConditions::InitConditions(string name, double f(vector<double>&, double, double), double T, vector<int> symmetry) :
	T(T), symmetry_vec(symmetry), init_cond_f(f), symmetry_nodes_num(symmetry[0] * symmetry[1] * symmetry[2]){
}

double InitConditions::f(vector<double>& vec) {
	return init_cond_f(vec, T, dz_cut);
}
vector<int> InitConditions::get_symmetry_vec() {
	return symmetry_vec;
}
vector<vector<int>> InitConditions::compute_symmetry_nodes(int N, int n_x, int n_y, int n_z) {
	vector<vector<int>> sym_int_nodes_vec(symmetry_nodes_num);
	int N_1 = symmetry_vec[0];
	int N_2 = symmetry_vec[1];
	int N_3 = symmetry_vec[2];
	vector<int> buff_1 = { n_x, N - n_x - 1 };
	vector<int> buff_2 = { n_y, N - n_y - 1 };
	vector<int> buff_3 = { n_z, N - n_z - 1 };

	for (int i = 0; i < sym_int_nodes_vec.size(); i++) {
		sym_int_nodes_vec[i] = { buff_1[i % N_1], buff_2[(i / N_1) % N_2], buff_3[(i / (N_1 * N_2)) % N_3] };
	}
	return sym_int_nodes_vec;

}
int InitConditions::get_symmetry_nodes_num() {
	return symmetry_nodes_num;
}
double InitConditions::get_temperature() {
	return T;
}
string InitConditions::get_f_name() {
	return name;
}
int InitConditions::symm_x(int N) {
	return N / symmetry_vec[0];
}
int InitConditions::symm_y(int N) {
	return N / symmetry_vec[1];
}
int InitConditions::symm_z(int N) {
	return N / symmetry_vec[2];
}
