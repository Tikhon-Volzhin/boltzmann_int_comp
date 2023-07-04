#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <string> 
#include <thread>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <iomanip>
#include <condition_variable>

#include <mutex>
#include <time.h>


#include "headers/VectorOverloading.h"
#include "headers/Net.h"
#include "headers/InitConditions.h"
#include "headers/AdditionalFunctions.h"

using namespace std;
int complited_threads;
int Num_threads = 10;

condition_variable cv;
Net::Net(int N, int num_steps, InitConditions & init_cond, double tau)
	:N(N), init_cond(init_cond), num_steps(num_steps), tau(tau), 
	dz_delta(2 * dz_cut / N), delt_dz3(dz_delta* dz_delta* dz_delta), 
	N_slope_dz(N/2 + 2){
	distribution_net = Vec_4d(num_steps + 1, Vec_3d(N, vector<vector<double>>(N, vector<double>(N))));;

	T = init_cond.get_temperature();
	
	Symm_x = init_cond.symm_x(N);
	Symm_y = init_cond.symm_y(N);
	Symm_z = init_cond.symm_z(N);
	
	

	mt19937 gen(rd());
	uniform_real_distribution<> dis(0.0, 1.0);
	this->gen = gen;
	this->dis = dis;

	string T_str = to_string(T);
	T_str.resize(3);
	net_file_prefix = "Log\\N" + to_string(N) + "_T" + T_str + "_" + init_cond.get_f_name();
	log_prefix = init_cond.get_f_name() + " N: " + to_string(N) + " num_steps: " + to_string(num_steps) + " T:" + T_str;

	create_log();
	create_velocity_array();
	gen_distribution();
	
	V_sph = N_0 * delt_dz3;
	V_sph_approx = (4.0 * pi * dz_cut * dz_cut * dz_cut / 3.0);

	num_integral_factor = (V_sph /sqrt(32) * N_0 * b_max * b_max) * tau;
	an_integral_factor = V_sph_approx / sqrt(2);     
}
void Net::gen_distribution() {
	double part_function = 0;
	N_0 = 0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				vector<int> ijk = { i,j,k };
				vector<double> dz_ijk = compute_dz(ijk);
				if ((dz_ijk ^ dz_ijk) < dz_cut_sq) N_0++;
				distribution_net[0][i][j][k] = init_cond.f(dz_ijk);
				if (distribution_net[0][i][j][k] < 0) cout << "jkhykj" << endl;
				part_function += distribution_net[0][i][j][k] * delt_dz3;
			}
		}
	}
	distribution_net[0] = distribution_net[0] / part_function;
	cout << "Net was generated" << endl;
}
void Net::create_velocity_array() {
	velocity_int_to_db = vector <double> (2*N + 4);

	for (int i = 0; i < N_slope_dz; i++) {
		velocity_int_to_db[i] = 2 * dz_cut;
	}

	for (int i = N_slope_dz; i < N_slope_dz + N; i++) {
		velocity_int_to_db[i] = dz_delta * ((double)(i - N/2 - 2) + 0.5) - dz_cut;
	}

	for (int i = N_slope_dz + N; i < 2*N + 4; i++) {
		velocity_int_to_db[i] = 2 * dz_cut;
	}
}
double Net::F(vector<int>& vec, int t) {
	return distribution_net[t][vec[0]][vec[1]][vec[2]];
}
double Net::F(vector<double>& buff, int t) {
	vector <int> vec = compute_int_coord(buff);
	return distribution_net[t][vec[0]][vec[1]][vec[2]];
}
void Net::F_update(vector<int>& vec, double delta_I, int t) {
	distribution_net[t + 1][vec[0]][vec[1]][vec[2]] += delta_I;
}


vector<double> Net::gen_rand_vec(int s) {
	vector<double> rand_vec(s);
	for (int i = 0; i < s; i++) {
		rand_vec[i] = dis(gen);
	}
	return rand_vec;
}
double Net::get_theta(double b_sq) {
	if (b_sq <= 1) {
		return 2 * acos(sqrt(b_sq));
	}
	return 0;
}

vector<vector<double>> Net::transform_velocities(vector<double>& dz_1, vector<double>& dz_2, vector<double>& angles) {
	vector<double> g = dz_2 - dz_1;
	double g_x = g[0], g_y = g[1], g_z = g[2], theta = angles[0], eps = angles[1];
	double g_norm = sqrt(g ^ g);

	double g_xy = sqrt(g_x * g_x + g_y * g_y);
	double cos_theta = cos(theta), cos_eps = cos(eps), sin_theta = sin(theta), sin_eps = sin(eps);
	
	if (!are_same_db_weak(g_xy, 0.0)) {
		g[0] = g_x * cos_theta - g_x * g_z / g_xy * cos_eps * sin_theta + g_norm * g_y / g_xy * sin_eps * sin_theta;
		g[1] = g_y * cos_theta - g_y * g_z / g_xy * cos_eps * sin_theta - g_norm * g_x / g_xy * sin_eps * sin_theta;
		g[2] = g_z * cos_theta + g_xy * cos_eps * sin_theta;
	}
	else {
		g[0] = g_norm * sin_eps * sin_theta;
		g[1] = g_norm * cos_eps * sin_theta;
		g[2] = g_norm * cos_theta;

	}
	vector<vector<double>> dz_new = { ((dz_1 + dz_2) / 2) - (g / 2), ((dz_1 + dz_2) / 2) + (g / 2) };
	return dz_new;

}

vector<double> Net::compute_dz(vector<int>& vec) {
	double i = velocity_int_to_db[vec[0] + N_slope_dz];
	double j = velocity_int_to_db[vec[1] + N_slope_dz];
	double k = velocity_int_to_db[vec[2] + N_slope_dz];
	return vector<double> { i, j, k };
}

vector<int> Net::compute_int_coord(vector<double>& vec) {
	return round_vec((vec / dz_delta) + (0.5 * N - 0.5) * ones_vec);
}
vector<int> Net::compute_int_center(vector<double>& vec) {
	return floor_vec((vec / dz_delta) + (0.5 * N - 0.5)* ones_vec);
}
vector<double> Net::approx_dz(vector<double>& vec) {
	vector <int> buff = compute_int_coord(vec);
	return compute_dz(buff);
}

void Net::find_nodes(vector <double>& dz_alpha, vector <double>& dz_beta, vector <int>& lambda, vector <int>& mu, vector <int>& lambda_add_s, vector <int>& mu_sub_s, double* r_v, bool* drop) {

	vector<int> alpha_near = compute_int_coord(dz_alpha);
	vector<int> alpha_center = compute_int_center(dz_alpha);

	vector<int> beta_near = compute_int_coord(dz_beta);
	
	if ((compute_dz(alpha_near) ^ compute_dz(alpha_near)) >= dz_cut_sq || (compute_dz(beta_near) ^ compute_dz(beta_near)) >= dz_cut_sq) {
		*drop = true;
		return;
	}

	vector<double> v_mass_center = (dz_alpha + dz_beta) / 2;
	vector<int> e_near = alpha_near - alpha_center;
	vector<vector<double>> dz_ijk(8, vector<double>(3));
	vector<double> E_ijk(8);

	double E_near, E_0 = ((dz_alpha - dz_beta) ^ (dz_alpha - dz_beta)) / 4.0;
	for (int i = 0; i < 8; i++) {
		vector<int> buff = alpha_center + n_ijk[i];
		dz_ijk[i] = compute_dz(buff);
		E_ijk[i] = (v_mass_center - dz_ijk[i]) ^ (v_mass_center - dz_ijk[i]);
		if (n_ijk[i] == e_near) {
			E_near = E_ijk[i];
		}
	}
	double min_dist_sq = dz_cut_sq;
	if (are_same_db_weak(E_0,E_near)) {
		lambda = alpha_near;
		lambda_add_s = alpha_near;
		mu = beta_near;
		mu_sub_s = beta_near;
		*r_v = 1.0;
	}
	else if (E_0 < E_near) {
		int ind = -1;
		double dist;
		
		for (int i = 0; i < 8; i++) {
			dist = (dz_ijk[i] - dz_alpha) ^ (dz_ijk[i] - dz_alpha);
			
			if (E_0 >= E_ijk[i] && min_dist_sq > dist && (dz_ijk[i] ^ dz_ijk[i]) < dz_cut_sq) {
			
				ind = i;
				min_dist_sq = dist;

			}
		}
		if (ind == -1) {
			*drop = true;
			return;
		}
		
		lambda_add_s = alpha_near;
		mu_sub_s = beta_near;
		lambda = alpha_center + n_ijk[ind];
		mu = (mu_sub_s - lambda) + lambda_add_s;
		if ((compute_dz(mu) ^ compute_dz(mu)) >= dz_cut_sq) {
			*drop = true;
			return;
		}

		double E_0v = E_0;
		double E_1v = E_ijk[ind];
		double E_2v = E_near;
		
		*r_v = (E_0v - E_1v) / (E_2v - E_1v);
	}
	else {

		int ind = -1;
		double dist;

		for (int i = 0; i < 8; i++) {
			dist = (dz_ijk[i] - dz_alpha) ^ (dz_ijk[i] - dz_alpha);
			if (E_0 <= E_ijk[i] && min_dist_sq > dist && (dz_ijk[i] ^ dz_ijk[i]) < dz_cut_sq) {
				ind = i;
				min_dist_sq = dist;
			}
		}
		if (ind == -1) {
			*drop = true;
			return;
		}
		lambda = alpha_near;
		mu = beta_near;
		lambda_add_s = alpha_center + n_ijk[ind];
		mu_sub_s = (lambda + mu) - lambda_add_s;
		if ((compute_dz(mu_sub_s) ^ compute_dz(mu_sub_s)) >= dz_cut_sq) {
			*drop = true;
			return;
		}
		double E_0v = E_0;
		double E_1v = E_near;
		double E_2v = E_ijk[ind];

		*r_v = (E_0v - E_1v) / (E_2v - E_1v);
		
	}
}
void Net::make_net_iter(int thread_ind,int t, vector<double>& rand_vec, pair<int,int> thread_split) {
	int thread_ind_start = thread_split.first;
	int thread_ind_end = thread_split.second;
	int num_collisions = thread_ind_end - thread_ind_start;
	
		vector<double> collision(8);

		vector<vector<double>> dz_transf_1;
		vector<vector<double>> dz_transf_2;
		vector<vector<double>> dz_clean_1;
		vector<vector<double>> dz_clean_2;
		vector<double> b_v;
		int approx_dropped_size = num_collisions;
		dz_transf_1.reserve(approx_dropped_size);
		dz_transf_1.reserve(approx_dropped_size);
		dz_clean_1.reserve(approx_dropped_size);
		dz_clean_2.reserve(approx_dropped_size);
		b_v.reserve(approx_dropped_size);
		for (int k = thread_ind_start; k < thread_ind_end; k++) {
			for (int s = 0; s < 8; s++) {
				double buff = (((long long)(k + 1) * num_corobov_params[s]) % num_corobov_prime) / (double)num_corobov_prime + rand_vec[s];
				collision[s] = buff - floor(buff);

				if (s < 6) {

					collision[s] = 2 * collision[s] * dz_cut - dz_cut;
				}
				else if (s == 7) {
					collision[s] = collision[s] * 2 * pi;
				}
			}

			vector<double> dz_1 = { collision[0] , collision[1] , collision[2] };
			vector<double> dz_2 = { collision[3] , collision[4] , collision[5] };
			dz_1 = approx_dz(dz_1);
			dz_2 = approx_dz(dz_2);


			if ((dz_1 ^ dz_1) <= dz_cut_sq && (dz_2 ^ dz_2) <= dz_cut_sq && (!(dz_1 == dz_2))) {
				dz_clean_1.push_back(dz_1);
				dz_clean_2.push_back(dz_2);
				vector<double> theta_eps_angles = { get_theta(collision[6]), collision[7] };
				vector<vector<double>> dz_new = transform_velocities(dz_1, dz_2, theta_eps_angles);
				dz_transf_1.push_back(dz_new[0]);
				dz_transf_2.push_back(dz_new[1]);

			}

		}
		int new_size = dz_clean_1.size();
		
		vector<int> shuffle_net(new_size);

		for (int i = 0; i < new_size; i++) {
			shuffle_net[i] = i;
		}

		shuffle(shuffle_net.begin(), shuffle_net.end(), gen);
		vector <int> lambda(3), mu(3), lambda_add_s(3), mu_sub_s(3);
		
		vector<int> alpha_v, beta_v;
		Mutex.lock();
		count_entire_collisions += new_size;
		Mutex.unlock();
		{
			std::unique_lock<std::mutex> lock(Mutex);
			complited_threads++;
		}
		cv.notify_all();
		{
			std::unique_lock<std::mutex> lock(Mutex);
			cv.wait(lock, []() { return (complited_threads==Num_threads); });
		}

		for (int sh_ind = 0; sh_ind < new_size; sh_ind++) {
			int i = shuffle_net[sh_ind];
			vector <double> dz_alpha_v = dz_transf_1[i], dz_beta_v = dz_transf_2[i];
			
			alpha_v = compute_int_coord(dz_clean_1[i]);
			beta_v = compute_int_coord(dz_clean_2[i]);
			
			double r_v;
			
			bool drop = false;
			find_nodes(dz_alpha_v, dz_beta_v, lambda, mu, lambda_add_s, mu_sub_s, &r_v, &drop);
			
			if (!drop) {
				double C_Omega_v;
				double alpha_beta_dist = sqrt((dz_alpha_v - dz_beta_v) ^ (dz_alpha_v - dz_beta_v));
				double C_dist = num_integral_factor * alpha_beta_dist/(double)count_entire_collisions;

				Mutex.lock();
				double f_lambda_mu = F(lambda, t + 1) * F(mu, t + 1);
				double f_lambda_mu_s = F(lambda_add_s, t + 1) * F(mu_sub_s, t + 1);
				double f_alpha_beta = F(alpha_v, t + 1) * F(beta_v, t + 1);
				
				//cout << f_lambda_mu << " " << f_lambda_mu_s << endl;;
				if (!are_same_db_weak(f_lambda_mu,0.0)) {
					C_Omega_v = C_dist * (f_lambda_mu * pow(f_lambda_mu_s / f_lambda_mu, r_v) - f_alpha_beta);
				}
				else {
					C_Omega_v = C_dist * (pow(f_lambda_mu, 1 - r_v) * pow(f_lambda_mu_s, r_v) - f_alpha_beta);
				}
				
				if ((F(alpha_v, t + 1) + C_Omega_v) < 0 ||
					(F(beta_v, t + 1) + C_Omega_v) < 0 ||
					(F(lambda, t + 1) + (r_v - 1.0) * C_Omega_v) < 0 ||
					(F(mu, t + 1) + (r_v - 1.0) * C_Omega_v) < 0 ||
					(F(lambda_add_s, t + 1.0) - r_v * C_Omega_v) < 0 ||
					(F(mu_sub_s, t + 1.0) - r_v * C_Omega_v) < 0) {
					continue;
				}
				F_update(alpha_v, C_Omega_v, t);
				F_update(beta_v, C_Omega_v, t);
				F_update(lambda, (r_v - 1.0) * C_Omega_v, t);
				F_update(mu, (r_v - 1.0) * C_Omega_v, t);
				F_update(lambda_add_s, -r_v * C_Omega_v, t);
				F_update(mu_sub_s, -r_v * C_Omega_v, t);
				count_coll++;
				if (count_coll % (int)(num_corobov_prime / 100) == 0) cout << "*";
				Mutex.unlock();

			}
		}
}

void Net::compute_num_I() {
	cout << "Numerical Integral computing" << endl;

	clock_t whole_calc = clock();

	vector<pair<int,int>> thread_split_arr(Num_threads + 1, pair<int, int>());
	int c_step = int(floor(num_corobov_prime / (double)Num_threads));
	for (int i = 1; i <= Num_threads; i++) {
		thread_split_arr[i].first = c_step * (i - 1);
		thread_split_arr[i].second = c_step * (i);
		if (i == Num_threads) thread_split_arr[i].second = num_corobov_prime;
	}

	for (int t = 0; t < num_steps; t++) {
		count_coll = 0;
		complited_threads = 0;
		count_entire_collisions = 0;

		clock_t calc = clock();
		distribution_net[t + 1] = distribution_net[t];
		vector<double> rand_vec = gen_rand_vec(8);

		vector<thread> thread_arr;
		for (int i = 1; i <= Num_threads; i++) {
			thread_arr.push_back(thread(&Net::make_net_iter,this,i, t, ref(rand_vec), thread_split_arr[i]));
		}
		
		for (int j = 0; j < Num_threads; j++) {
			thread_arr[j].join();
		}
		cout << endl << "Step: " << t + 1 << ". " << count_coll << " collisions" << endl;

		printf("Time taken: %.2fs\n", (double)(clock() - calc) / CLOCKS_PER_SEC);

		double I_norm = (distribution_net[t + 1] - distribution_net[0]) ^ (distribution_net[t + 1] - distribution_net[0]);
		cout << "Integral norm: " << sqrt( I_norm * delt_dz3 ) / tau / (t + 1) << endl;
		cout << "Entire collisions: " << count_entire_collisions << endl;
		
	}
	printf("####Overall time taken: %.2fs####\n\n\n", (double)(clock() - whole_calc) / CLOCKS_PER_SEC);

	double time = tau * num_steps;
	numerical_Integral = (distribution_net[num_steps] - distribution_net[0]) / time;
	
	compute_sym_numerical_Integral();

	string f_name = net_file_prefix +"_" + to_string(num_steps)+ "_Num_vec.txt";
	write_3dvec_in_file(numerical_Integral, f_name);
	is_num_integral_created = true;
	int t = 0;
	for (auto a : distribution_net) {
		string name = "Dist" + to_string(N)+"//" + to_string(t) + ".txt";
		write_3dvec_in_file(a, name);
		t++;
	}
	test_num_integral();
}
void Net::analytical_thread_function(pair<int,int> thread_split) {
	for (int i = thread_split.first; i < thread_split.second; i++) {
		for (int j = 0; j < init_cond.symm_y(N); j++) {
			for (int k = 0; k < init_cond.symm_z(N); k++) {
				double I_ijk = 0;
				vector<int> ijk = { i,j,k };
				vector<double> dz_node = compute_dz(ijk);
				double buff = (dz_node ^ dz_node);
				if (buff >= dz_cut_sq) {
					continue;
				}
				int count_collisions_an = 0;
				vector<double> collision(5);
				for (int p = 0; p < an_corobov_prime; p++) {
					
					for (int s = 0; s < 5; s++) {
						double buff = (((long long)(p + 1) * an_corobov_params[s]) % an_corobov_prime) / (double)an_corobov_prime;
						collision[s] = buff - floor(buff);

						if (s < 3) {

							collision[s] = 2 * collision[s] * dz_cut - dz_cut;
						}
						else if (s == 4) {
							collision[s] = collision[s] * 2 * pi;
						}
					}

					vector<double> dz = { collision[0] , collision[1] , collision[2] };
					if ((dz ^ dz) >= dz_cut_sq || (dz_node == dz))
						continue;
					count_collisions_an++;
					vector<double> theta_eps_angles = { get_theta(collision[3]), collision[4] };
					vector<vector<double>> dz_new = transform_velocities(dz_node, dz, theta_eps_angles);
					if ((dz_new[0] ^ dz_new[0]) >= dz_cut_sq || (dz_new[1] ^ dz_new[1]) >= dz_cut_sq) {
						continue;
					}

					double f_1_trans = init_cond.f(dz_new[0]);
					double f_2_trans = init_cond.f(dz_new[1]);
					double f_init = init_cond.f(dz_node);
					double f_node = init_cond.f(dz);

					double dist = sqrt((dz_node - dz) ^ (dz_node - dz));

					I_ijk += (f_1_trans * f_2_trans - f_init * f_node) * dist;
				}
				Mutex.lock();
				analytical_Integral[i][j][k] = I_ijk / (double)count_collisions_an;
				Mutex.unlock();
			}
		}
		Mutex.lock();
		cout << "#";
		Mutex.unlock();
	}

}
void Net::compute_an_I() {
	
	int analitycal_num_threads = 10;

	if (Symm_x < 10) {
		analitycal_num_threads = Symm_x;
	}

	vector<pair<int, int>> thread_split_arr(analitycal_num_threads + 1, pair<int,int>());
	int c_step = Symm_x / analitycal_num_threads;
	for (int i = 1; i <= analitycal_num_threads; i++) {
		thread_split_arr[i].first = c_step * (i - 1);
		thread_split_arr[i].second = c_step * (i);
		if (i == analitycal_num_threads) thread_split_arr[i].second = Symm_x;
	}

	analytical_Integral = Vec_3d(Symm_x, vector<vector<double>>(Symm_y, vector<double>(Symm_z, 0)));

	cout << "Analytical Integral computing" << endl;
	clock_t whole_calc = clock();

	vector<thread> thread_arr;
	for (int i = 1; i <= analitycal_num_threads; i++) {
		thread_arr.push_back(thread(&Net::analytical_thread_function, this, thread_split_arr[i]));
	}

	for (int j = 0; j < analitycal_num_threads; j++) {
		thread_arr[j].join();
	}

	printf("####Overall time taken: %.2fs####\n\n", (double)(clock() - whole_calc) / CLOCKS_PER_SEC);
	analytical_Integral = an_integral_factor * analytical_Integral;
	ofstream I_f_cons("Log\\Log.txt", ios::app);
	if (I_f_cons.is_open()) {
		double I_val = 0;
		vector<double> I_momentum = { 0.0, 0.0, 0.0 };
		double I_energy = 0;
		double I_sq_norm = 0;

		for (int i = 0; i < Symm_x; i++) {
			for (int j = 0; j < Symm_y; j++) {
				for (int k = 0; k < Symm_z; k++) {
					vector<vector<int>> symm_nodes = init_cond.compute_symmetry_nodes(N, i, j, k);
					for (int l = 0; l < symm_nodes.size(); l++) {
						vector<double> momentum = compute_dz(symm_nodes[l]);
						I_sq_norm += analytical_Integral[i][j][k] * analytical_Integral[i][j][k] * delt_dz3;
						I_val += analytical_Integral[i][j][k] * delt_dz3;
						I_momentum = I_momentum + (analytical_Integral[i][j][k] * delt_dz3) * momentum;
						I_energy += (analytical_Integral[i][j][k] * delt_dz3) * (momentum ^ momentum);
					}
				}
			}
		}

		double I_norm = sqrt(I_sq_norm);
		I_val /= I_norm;
		I_momentum = I_momentum / I_norm;
		I_energy /= I_norm;

		I_f_cons << endl << log_prefix << endl;
		I_f_cons << "Analytical: " << I_norm << " concentration: " << I_val << " momentum : " << I_momentum << " energy: " << I_energy << endl;
		cout << "Analytical: " << I_norm << " concentration: " << I_val << " momentum : " << I_momentum << " energy: " << I_energy << endl;
	}
	I_f_cons.close();
	string f_name = net_file_prefix + "_An_vec.txt";
	write_3dvec_in_file(analytical_Integral, f_name);

	is_an_integral_created = true;
}

void Net::compute_sym_numerical_Integral() {
	sym_numerical_Integral = Vec_3d(Symm_x, vector<vector<double>>(Symm_y, vector<double>(Symm_z, 0)));
	for (int i = 0; i < Symm_x; i++) {
		for (int j = 0; j < Symm_y; j++) {
			for (int k = 0; k < Symm_z; k++) {
				vector<vector<int>> sym_nodes = init_cond.compute_symmetry_nodes(N, i, j, k);
				for (int l = 0; l < sym_nodes.size(); l++) {
					int _i = sym_nodes[l][0], _j = sym_nodes[l][1], _k = sym_nodes[l][2];
					sym_numerical_Integral[i][j][k] += numerical_Integral[_i][_j][_k];
				}
			}
		}
	}
	sym_numerical_Integral = sym_numerical_Integral/init_cond.get_symmetry_nodes_num();
}
void Net::compute_convergence_delta() {

	if (!is_an_integral_created)
		read_an_I();
	if (!is_num_integral_created) {
		read_num_I();
		compute_sym_numerical_Integral();
	}
	double an_norm = sqrt((analytical_Integral ^ analytical_Integral) * delt_dz3);
	double num_norm = sqrt((sym_numerical_Integral ^ sym_numerical_Integral) * delt_dz3);
	cout << endl << endl << "############ TEST 3 #############" << endl;
	cout << "##### Convergence Delta #####" << endl;
	cout << "an_norm: " <<an_norm << " num_norm: " << num_norm << " An/Num: " << an_norm/ num_norm << endl;
	cout << "mean: " << mean_3d_vec(sym_numerical_Integral/analytical_Integral ) << " std: " << std_3d_vec(sym_numerical_Integral / analytical_Integral) << endl;
	double convergence_delta = sqrt(((analytical_Integral - sym_numerical_Integral) ^ (analytical_Integral - sym_numerical_Integral)) * delt_dz3)/an_norm;
	double error = sqrt(((analytical_Integral - sym_numerical_Integral) ^ (analytical_Integral - sym_numerical_Integral)) * delt_dz3);
	cout << "convergence delta: " << convergence_delta << " error: " << error <<endl;

	ofstream file("Log\\Log.txt", ios::app);
	if (file.is_open()) {
		file << endl << log_prefix << endl;
		file << "convergence delta: " << convergence_delta << endl;
	}
	file.close();
}

void Net::read_an_I() {
	string f_name = net_file_prefix + "_An_vec.txt";
	analytical_Integral = read_3dvec_from_file(f_name);
}
void Net::read_num_I() {
	string f_name = net_file_prefix + "_" + to_string(num_steps) + "_Num_vec.txt";
	numerical_Integral = read_3dvec_from_file(f_name);
}

void Net::test_num_integral() {
	cout << "############ TEST 1 #############" << endl;
	cout << "##### Conserving quantities #####" << endl;
	string start_end_prefix = "Start";
	for (int t = 0; t < num_steps + 1; t = t + num_steps) {
		double part_function = 0;
		vector<double> total_momentum = { 0.0,0.0,0.0 };
		double total_energy = 0;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				for (int k = 0; k < N; k++) {
					part_function += distribution_net[t][i][j][k] * delt_dz3;
					vector<int> ijk = { i,j,k };
					total_momentum = total_momentum + (distribution_net[t][i][j][k] * delt_dz3 * compute_dz(ijk));
					total_energy = total_energy + (distribution_net[t][i][j][k] * delt_dz3 * (compute_dz(ijk) ^ compute_dz(ijk)));
				}
			}
		}
		cout << start_end_prefix << " Distribution: " << part_function << " Total momentum: [" << total_momentum[0] << "," << total_momentum[1] << "," << total_momentum[2] << "] Energy: " << total_energy << endl;
		start_end_prefix = "End";
	}
	cout << "############ TEST 2 #############" << endl;
	cout << "######## Mirror symmetry ########" << endl;

	double I_x = 0.0, I_y = 0.0, I_z = 0.0;
	double I_abs = (numerical_Integral ^ numerical_Integral);
	cout << "Numerical Integral = " << sqrt(I_abs * delt_dz3) << endl;

	for (int j = 0; j < N; j++) {
		for (int k = 0; k < N; k++) {
			for (int i = 0; i < N / 2; i++) {
				double diff = numerical_Integral[i][j][k] - numerical_Integral[N - 1 - i][j][k];
				I_x += diff * diff;
			}
		}
	}

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N / 2; j++) {
			for (int k = 0; k < N; k++) {
				double diff = numerical_Integral[i][j][k] - numerical_Integral[i][N - 1 - j][k];
				I_y += diff * diff;
			}
		}
	}

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N / 2; k++) {
				double diff = numerical_Integral[i][j][k] - numerical_Integral[i][j][N - 1 - k];
				I_z += diff * diff;
			}
		}
	}



	cout << "dI_x = " << sqrt(2 * I_x / I_abs) << "; dI_y = " << sqrt(2 * I_y / I_abs) << "; dI_z = " << sqrt(2 * I_z / I_abs) << endl;

	ofstream file("Log\\Log.txt", ios::app);
	if (file.is_open()) {
		file << endl << log_prefix << endl;
		file << "Numerical Integral = " << sqrt(I_abs * delt_dz3) << " dI_x = " << sqrt(2 * I_x / I_abs) << "; dI_y = " << sqrt(2 * I_y / I_abs) << "; dI_z = " << sqrt(2 * I_z / I_abs) << endl;
	}
	file.close();
}
void Net::create_log() {
	ifstream file("Log\\Log.txt", ios::app);
	if (!file.is_open()) {
		ofstream file_log("Log\\Log.txt");
		file_log.close();
	}
	file.close();
}