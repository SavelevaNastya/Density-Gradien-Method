#pragma once
#include <math.h>
#include <vector>
#include <map>
#include <mpi.h>
#include <fstream>
#include <sstream>
class solver {
private:
	std::map<std::string, std::vector<double>> settings; //N, M, delta_t, delta_z, eps, Temp, b, omega, Tc[], Pc[]

	std::vector<double> Tr;
	std::vector<double> alpha_Tr_omega;
	std::vector<double> b_i;
	std::vector<double> a_i;
	std::vector<std::vector<double>> c, Mcoeff, Acoeff;

	int numOfProc;
	int sendcount;

	double delta_1 = 1 - sqrt(2);
	double delta_2 = 1 + sqrt(2);
	double R = 8.31446261815324;
	double m_omega;

	std::vector<std::vector<double>> n; // first index - number of comp.; second - spatial step
	std::vector<double> f, x; std::vector<double> A_value, Arecv_value;; // SLAE A*x = f
	
	std::vector<int> A_strIdx, Arecv_strIdx;
	std::vector<int> A_colIdx, Arecv_colIdx;

	void subtract_vec();
	void abs_subtract_vec();
	void update_X();
	void devide();
	void assignment();
	double normaInf();
	void mult_matvec(int rank);
	std::vector<double> TempX;
	std::vector<double> Ax;
	std::vector<double> rezX;

	double df_bulk(int, int);
	void RightSideVector();
	void matrixA();
	void Mcoeff_init();
	void Acoeff_init();
	void Ccoeff_init();
	void n_initDist_init();
public:
	solver() {};
	void initialization(int rank);
	void Jacobi(int rank);
	void printA();
	void printAnswer(int rank);
};

