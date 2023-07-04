#ifndef ADDITIONAL_FUNCTIONS
#define ADDITIONAL_FUNCTIONS
#define Vec_4d vector<vector<vector<vector<double>>>>
#define Vec_3d vector<vector<vector<double>>>
using namespace std;


bool are_same_db_strict(double a, double b);
bool are_same_db_weak(double a, double b);
bool are_same_db_e5(double a, double b);


vector<int> round_vec(const vector<double>& vec);
vector<int> floor_vec(const vector<double>& vec); 

ostream& operator<< (ostream& out, const vector<int>& vec);
ostream& operator<< (ostream& out, const vector<double>& vec);

void write_3dvec_in_file(const Vec_3d& vec, string name, vector<int> cut = { 1,1,1 });
Vec_3d read_3dvec_from_file(string name);

double mean_3d_vec(const Vec_3d& vec);
double std_3d_vec(const Vec_3d& vec);



#endif