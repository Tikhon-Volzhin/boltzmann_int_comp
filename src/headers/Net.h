#ifndef NET
#define NET
#include "headers/InitConditions.h"
#define Vec_4d vector<vector<vector<vector<double>>>>
#define Vec_3d vector<vector<vector<double>>>
using namespace std;
class Net
{
public:
	Net(int N, int num_steps, InitConditions & init_cond, double tau = pow(10, -6));

	void gen_distribution();
	void compute_num_I();
	void compute_an_I();
	void compute_convergence_delta();
	void test_num_integral();
	void read_an_I();
	void read_num_I();
	
	
private:
	double F(vector<int>& vec, int t);
	double F(vector<double>& buff, int t);
	void F_update(vector<int>& vec, double delta_I, int t);
	vector<double> gen_rand_vec(int s);
	double get_theta(double b_sq);
	vector<vector<double>> transform_velocities(vector<double>& dz_1, vector<double>& dz_2, vector<double>& angles);
	void find_nodes(vector <double>& dz_alpha, vector <double>& dz_beta, vector <int>& lambda, vector <int>& mu, vector <int>& lambda_add_s, vector <int>& mu_sub_s, double* r_v, bool* drop);
	void create_velocity_array();
	void make_net_iter(int thread_ind,int t, vector<double>& rand_vec, pair<int, int> thread_split);
	void compute_sym_numerical_Integral();
	vector<double> compute_dz(vector<int>& vec);
	vector<int> compute_int_coord(vector <double>& vec);
	vector<int> compute_int_center(vector<double>& vec);
	vector<double> approx_dz(vector<double>& vec);
	void analytical_thread_function(pair<int, int> thread_split);
	void create_log();

	string net_file_prefix;
	string log_prefix;
	InitConditions init_cond;

	int num_corobov_prime = 4000037;
	vector<int> num_corobov_params = { 1,72362, 210611, 92212, 583028, 681897, 2974319, 1680656 };

	int an_corobov_prime = 200003;
	vector<int> an_corobov_params = { 1, 58186, 159815, 56108, 51119 };

	int N, num_steps;

	double b_max = 1.0, tau;
	double dz_cut = 4.8, dz_cut_sq = dz_cut* dz_cut;
	double dz_delta, delt_dz3;
	double num_integral_factor;
	double an_integral_factor;
	double T;
	
	int Symm_x, Symm_y, Symm_z;
	
	bool is_num_integral_created = false;
	bool is_an_integral_created = false;
	double V_sph, V_sph_approx;
	int N_0;

	Vec_4d distribution_net;
	Vec_3d analytical_Integral;
	Vec_3d numerical_Integral;
	Vec_3d sym_numerical_Integral;

	vector<double> velocity_int_to_db;
	vector<int> ones_vec = { 1,1,1 };
	vector<vector<int>> n_ijk = { { 0,0,0}, { 0,0,1 }, { 0,1,0 }, { 1,0,0 }, { 1,1,0 }, { 0,1,1 }, { 1,0,1 }, { 1,1,1 } };

	random_device rd;
	mt19937 gen;
	uniform_real_distribution<> dis;
	mutex Mutex;
	double pi = 3.14159265358979323846;
	int count_coll;
	int count_entire_collisions;
	int N_slope_dz;
};
#endif