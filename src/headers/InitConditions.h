#ifndef INIT_CONDITIONS
#define INIT_CONDITIONS
using namespace std;

class InitConditions {
public:
	InitConditions(string name,double f(vector<double>&, double, double), double T, vector<int> symmetry = { 1,1,1 });

	double f(vector<double>&);

	vector<int> get_symmetry_vec();
	vector<vector<int>> compute_symmetry_nodes(int N, int i, int j, int k);
	int get_symmetry_nodes_num();
	double get_temperature();
	string get_f_name();

	int symm_x(int N);
	int symm_y(int N);
	int symm_z(int N);

private:
	double T;
	double dz_cut = 4.8;
	string name;
	int symmetry_nodes_num;
	vector<int> symmetry_vec;
	double (*init_cond_f)(vector<double>& vec, double T, double dz_cut);
};
#endif