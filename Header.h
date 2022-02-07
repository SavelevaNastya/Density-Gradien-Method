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
	int RowNum;

	double delta_1 = 1 - sqrt(2);
	double delta_2 = 1 + sqrt(2);
	double R = 8.31446261815324;
	double m_omega;

	std::vector<std::vector<double>> n, mu; // first index - number of comp.; second - spatial step
	std::vector<double> f, x; std::map<std::pair<int, int>, double> A; // SLAE A*x = f
	
	void subtract_vec(std::vector<double> x, std::vector<double> y, std::vector<double>& rez, int rank);
	void abs_subtract_vec(std::vector<double> x, std::vector<double> y, std::vector<double>& z, int rank);
	void update_X(std::vector<double>& Xnew, std::vector<double>& X, int rank);
	void devide(std::vector<double>& X, std::map<std::pair<int, int>, double> A, int rank);
	void assignment(std::vector<double>& X, std::vector<double> F, int rank);
	double normaInf(std::vector<double> X, int rank);
	void mult_matvec_(std::map<std::pair<int, int>, double> A, std::vector<double> b, std::vector<double>& Ab, int rank);
	void mult_matvec(std::map<std::pair<int, int>, double> A, std::vector<double> b, std::vector<double>& Ab, int rank);
	double df_bulk(int, int);
	void RightSideVector();
	void matrixA();
public:
	solver() {};
	void initialization(int);
	void Jacobi(int rank);
};

