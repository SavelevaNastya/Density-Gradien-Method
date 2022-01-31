#pragma once
#include <math.h>
#include <vector>


class solver {
private:
	int M; //number of comp;
	int T; //time step
	int N; //spatial step
	double eps = 1e-7;

	double delta_1 = 1 - sqrt(2);
	double delta_2 = 1 + sqrt(2);
	double omega;
	double m_omega;
	void define_m_omega() {
		if (omega <= 0.491) {
			m_omega = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
		}
		else {
			m_omega = 0.379642 + 1.48503 * omega - 0.164423 * omega * omega + 0.016666 * omega * omega * omega;
		}
	};
	double R = 8.31446261815324; // Дж/(моль∙К).
	double Temp;
	double b;
	std::vector<double> Tc;
	std::vector<double> Pc;
	std::vector<double> Tr = Temp / Tc[i];
	std::vector<double> alpha_Tr_omega = (1 + m_omega * (1 - sqrt(Tr[i]))) * (1 + m_omega * (1 - sqrt(Tr[i])));
	std::vector<double> b_i = 0.0778 * R * Tc[i] / Pc[i];
	std::vector<double> a_i = 0.457235 * alpha_Tr_omega[i] * R * R * Tc[i] * Tc[i] / Pc[i];

	std::vector<std::vector<double>> n, mu; // first index - number of comp.; second - spatial step
	std::vector<std::vector<double>> c, Mcoeff, Acoeff;
	double D = 1.0, delta_z = D / (N - 1), delta_t = 1.0 / (T - 1);
	
	std::vector<double> f, x, x_new; std::vector<std::vector<double>> A; // SLAE A*x = f

	void subtract_vec(std::vector<double> x, std::vector<double> y, std::vector<double>& rez, int rank, int RowNum);
	void abs_subtract_vec(std::vector<double> x, std::vector<double> y, std::vector<double>& z, int rank, int RowNum);
	void update_X(std::vector<double>& Xnew, std::vector<double>& X, int rank, int RowNum);
	void devide(std::vector<double>& X, std::vector<std::vector<double>> A, int rank, int RowNum);
	void assignment(std::vector<double>& X, std::vector<double> F, int rank, int RowNum);
	double normaInf(std::vector<double> X, int rank, int RowNum);
	void mult_matvec_(std::vector<std::vector<double>> A, std::vector<double> b, std::vector<double>& Ab, int rank, int RowNum);
	void mult_matvec(std::vector<std::vector<double>> A, std::vector<double> b, std::vector<double>& Ab, int rank, int RowNum);
public:
	solver() {};
	void RightSideVector();
	void matrixA();
	double df_bulk(int, int);
    void Jacobi(std::vector<std::vector<double>> A, std::vector<double> F, std::vector<double>& X, int rank, int RowNum);
	
};

