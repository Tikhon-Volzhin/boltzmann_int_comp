#ifndef VEC_OVER
#define VEC_OVER

#define Vec_3d std::vector<vector<vector<double>>>

using namespace std;

vector<double> operator + (const vector<double>& v1, const vector<double>& v2);
vector<int> operator + (const vector<int>& v1, const vector<int>& v2);


vector<double> operator - (const vector<double>& v1, const vector<double>& v2);
vector<int> operator - (const vector<int>& v1, const vector<int>& v2);


double operator ^ (const vector<double>& v1, const vector<double>& v2);
int operator ^ (const vector<int>& v1, const vector<int>& v2);

vector<double> operator * (const vector<double>& v1, const vector<double>& v2);
vector<int> operator * (const vector<int>& v1, const vector<int>& v2);

vector<double> operator * (const double real_number, const vector<double>& vec);
vector<double> operator * (const double real_number, const vector<int>& vec);
vector<int> operator * (const int real_number, const vector<int>& vec);

vector<double> operator / (const vector<double>& vec, const double real_number);
vector<double> operator / (const vector<int>& vec, const double real_number);

bool operator == (const vector<int>& v1, const vector<int>& v2);
bool operator == (const vector<double>& v1, const vector<double>& v2);


Vec_3d operator + (const Vec_3d& v1, const Vec_3d& v2);
Vec_3d operator - (const Vec_3d& v1, const Vec_3d& v2);
Vec_3d operator / (const Vec_3d& v1, const Vec_3d& v2);


Vec_3d operator / (const Vec_3d& v1, double scalar);
Vec_3d operator / (const Vec_3d& v1, int scalar);

Vec_3d operator * (double scalar, const Vec_3d& v1);

double operator ^ (const Vec_3d& v1, const Vec_3d& v2);



#endif