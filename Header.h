#pragma once
#include <math.h>
#include <vector>
#define M 10 //number of comp
#define T 10 //time step
#define N 10 //spatial step
#define eps 1e-7 

class solver {
private:
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
	double Tc;
	double Pc;
	double R = 8.31446261815324; // Дж/(моль∙К).
	double Temp;
	double Tr = Temp / Tc;
	double alpha_Tr_omega = (1 + m_omega * (1 - sqrt(Tr))) * (1 + m_omega * (1 - sqrt(Tr)));
	double b = 0.0778 * R * Tc / Pc;
	double a = 0.457235 * alpha_Tr_omega * R * R * Tc * Tc / Pc;
	double b_i[M];  // ???????????????

	double n[M][T][N], n_new[M][T][N]; // first index - number of comp; second - time step; third - spatial step
	double c[M][M], Mcoeff[M][M], Acoeff[M][M];
	double D = 1.0, delta_z = D / (N - 1), delta_t = 1.0 / (T - 1);
	double W[M][M], W_inv[M][M], F[M];
public:
	solver() {};
	double f(int,int,int);
	double df_bulk(int, int, int);
	void init_cond();
	void WmatrixFormation(int, int, int);
	void matInverse(int, int, int);
	double* delta(int, int, int);
	void newtonMethod(int);
	double norm(int, int);
	void solve();
	void FvectorFormation(int, int);
};