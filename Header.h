#pragma once
#include <math.h>
#include <vector>
#include <map>
#include <mpi.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <F:/Eigen/Dense>
#include <F:/Eigen/Sparse>

using namespace Eigen;

class solver {
private:
	std::map<std::string, std::vector<double>> settings; //M, delta_t, delta_z, Temp, omega[], Tc[], Pc[]

	std::vector<double> Tr, alpha_Tr_omega, a_i, b_i, z, m_omega;
	std::vector<std::vector<double>> c, Mcoeff, Acoeff, Kcoeff, temp;

	double delta_1 = 1 - sqrt(2);
	double delta_2 = 1 + sqrt(2);
	double R = 8.31446261815324;

	std::vector<std::vector<double>> n, n_new; // first index - number of comp.; second - spatial step

	Eigen::VectorXd x, f;
	Eigen::SparseMatrix<double> A_value; //SLAE A* x = f
	double delta_z, delta_t, time;

	double df_bulk(int, int);
	void RightSideVector();
	void matrixA();
	void Mcoeff_init();
	void Kcoeff_init();
	void Lcoeff_init();
	void n_initDist_init();
	void fill_matrix_from_file(std::string path, char Matrix);
	void writeAnswer();
	void new_time_step();
	double f_sigma(int);
public:
	void initialization();
	void solve();
	void sigma();
};