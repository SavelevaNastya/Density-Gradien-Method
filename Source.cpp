#include "Header.h"
#include <iostream>

double solver::df_bulk(int compIdx, int spatialIdx) {
    define_m_omega();
    double n_sum = 0;
    double sum1 = 0, sum2 = 0;
    for (int i = 0; i < M; i++) {
        n_sum = n_sum + n[i][spatialIdx];
        sum1 += n[i][spatialIdx] * Acoeff[compIdx][i];
        sum2 += n[i][spatialIdx] * Acoeff[i][compIdx];
    }
    
    double rez = R * Temp * (b_i[compIdx] * n_sum / (1 - b * n_sum) + log(n[compIdx][spatialIdx]) - log(1 - b * n_sum))
        + ((sum1 + sum2)/(b*n_sum) - a_i[compIdx]*b_i[compIdx]/(b*b)) * log((b * delta_1 * n_sum + 1) / (b * delta_2 * n_sum + 1)) / (delta_2 - delta_1)
        - a_i[compIdx] * b_i[compIdx] * n_sum / (b * (delta_1 * b * n_sum + 1) * (delta_2 * b * n_sum + 1));
    return rez;
}

void solver::RightSideVector() {
    
    for (size_t i = 0; i < n.size(); i++) {
        for (size_t j = 0; j < n[i].size(); j++) {
            f.push_back(n[i][j]);
        }
    }

    for (size_t i = 0; i < mu.size(); i++) {
        for (size_t j = 0; j < mu[i].size(); j++) {
            f.push_back(df_bulk(i,j));
        }
    }
}

void solver::matrixA() {
    // n.size() = M, n[1].size() = N
    int M = n.size(), N = n[1].size();
    for (size_t k = 0; k < M; k++) { //номер блока по горизонтали 
        for (size_t p = 0; p < M; p++) { // номер блока по вертикали 
            for (size_t i = k * (N); i < (1 + k) * N; i++) { // номер строки в этом блоке
                for (size_t j = p * (N); j < (1 + p) * N; j++) { // номер столбца в этом блоке

                    if ( (i == k*N) || (i == (k+1)*N - 1) ) { 
                        A[i][j] = 0.0;
                        continue;
                    }
                    if ( (i - k * N) == (j - p * N) ) {
                        A[i][j] = 1.0;
                    }
                }
            }
        }
    }

    for (size_t k = 0; k < M; k++) { //номер блока по горизонтали 
        for (size_t p = 0; p < M; p++) { // номер блока по вертикали
            for (size_t i = k * N; i < (1 + k) * N; i++) { // номер строки в этом блоке
                for (size_t j = N * M + p * N; j < N * M + (1 + p) * N; j++) { // номер столбца в этом блоке
                    double coeff_1 = Mcoeff[k][p] * delta_t / (delta_z * delta_z);
                    if ( (i == k * N) || (i == (k + 1) * N - 1) ) { 
                        A[i][j] = 0.0;
                        continue;
                    }
                    if ( (i - k * N) == (j - N * M - p * N - 1) ) {
                        A[i][j] = coeff_1;
                    }
                    if ( (i - k * N) == (j - N * M - p * N) ) {
                        A[i][j] = -2*coeff_1;
                    }
                    if ((i - k * N) == (j - N * M - p * N + 1)) {
                        A[i][j] = coeff_1;
                    }
                }
            }
        }
    }
    
    for (size_t k = 0; k < M; k++) { //номер блока по горизонтали 
        for (size_t p = 0; p < M; p++) { // номер блока по вертикали 
            for (size_t i = N * M + k * N; i < N * M + (1 + k) * N; i++) { // номер строки в этом блоке
                for (size_t j = p * (N); j < (1 + p) * N; j++) { // номер столбца в этом блоке
                    double coeff_1 = c[k][p] / (delta_z * delta_z);
                    if ((i == N * M + k * N) || (i == (k + 1) * N - 1 + N * M)) {
                        A[i][j] = 0.0;
                        continue;
                    }
                    if ((i - k * N - N * M) == (j  - p * N - 1)) {
                        A[i][j] = coeff_1;
                    }
                    if ((i - k * N - N * M) == (j - p * N)) {
                        A[i][j] = -2 * coeff_1;
                    }
                    if ((i - k * N - N * M) == (j - p * N + 1)) {
                        A[i][j] = coeff_1;
                    }
                }
            }
        }
    }

    for (size_t k = 0; k < M; k++) { //номер блока по горизонтали 
        for (size_t p = 0; p < M; p++) { // номер блока по вертикали 
            for (size_t i = N * M + k * N; i < N * M + (1 + k) * N; i++) { // номер строки в этом блоке
                for (size_t j = N * M + p * N; j < N * M + (1 + p) * N; j++) { // номер столбца в этом блоке
                    
                    if ((i == N * M + k * N) || (i == (k + 1) * N - 1 + N * M)) {
                        A[i][j] = 0.0;
                        continue;
                    }
                    
                    if ((i - k * N - N * M) == (j - p * N - N * M)) {
                        A[i][j] = 1.0;
                    }else{ A[i][j] = 0.0; }
                    
                }
            }
        }
    }
}

double solver::normaInf(std::vector<double> X, int rank, int RowNum) {  
    //std::vector<double> abs_vec;
    //for (int compIdx = 0; compIdx < M; compIdx++) {
    //    //abs_vec[compIdx] = n[compIdx][spatialIdx] - n_new[compIdx][spatialIdx];
    // }
    //double rez = fabs(abs_vec[0]);
    //for (int i = 1; i < M; i++) {
    //    if (fabs(abs_vec[i]) > rez) {
    //        rez = fabs(abs_vec[i]);
    //    }
    //}
    //return rez;

    double norm = 0.0;
    for (int i = rank * RowNum; i < (rank + 1) * RowNum; i++)
    {
        MPI_Allreduce(&x[i], &norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    return norm;
}

void solver::mult_matvec(std::vector<std::vector<double>> A, std::vector<double> b, std::vector<double>& Ab, int rank, int RowNum) {

    for (int i = rank * RowNum; i < (rank + 1) * RowNum; i++) {
        for (int j = 0; j < N; j++) {
            if (i != j)
                Ab[i] += A[i][j] * b[j];
        }
    }
}

void solver::mult_matvec_(std::vector<std::vector<double>> A, std::vector<double> b, std::vector<double>& Ab, int rank, int RowNum) {

    for (int i = rank * RowNum; i < (rank + 1) * RowNum; i++) {
        for (int j = 0; j < N; j++) {
            Ab[i] += A[i][j] * b[j];
        }
    }
}

void solver::Jacobi(std::vector<std::vector<double>> A, std::vector<double> F, std::vector<double>& X, int rank, int RowNum) {
	
    std::vector<double> TempX;
    std::vector<double> Ax;
    std::vector<double> rezX;
	
	double norm = 0.0;
	int iter = 0;

	for (int i = 0; i < N; i++) {
		Ax[i] = 0.0;
		TempX[i] = 0.0;
		rezX[i] = 0.0;
	}

	do {
		iter++;
		assignment(TempX, F, rank, RowNum);
		mult_matvec(A, X, Ax, rank, RowNum);
		subtract_vec(TempX, Ax, TempX, rank, RowNum);
		devide(TempX, A, rank, RowNum);

		abs_subtract_vec(X, TempX, rezX, rank, RowNum);
		norm = normaInf(rezX, rank, RowNum);
		if (rank == 0) { std::cout << "NORMA = " << norm << " " << "ITER " << iter << std::endl; }

		update_X(TempX, X, rank, RowNum);

	} while (norm > eps);
	if (rank == 0) { std::cout << "	Number of iteration = " << iter << std::endl; }
}

void solver::subtract_vec(std::vector<double> x, std::vector<double> y, std::vector<double>& rez, int rank, int RowNum) {
    for (int i = rank * RowNum; i < (rank + 1) * RowNum; i++) {
        rez[i] = x[i] - y[i];
    }
}

void solver::abs_subtract_vec(std::vector<double> x, std::vector<double> y, std::vector<double>& z, int rank, int RowNum) {
    for (int i = rank * RowNum; i < (rank + 1) * RowNum; i++) {
        z[i] = fabs(x[i] - y[i]);
    }
}

void solver::update_X(std::vector<double>& Xnew, std::vector<double>& X, int rank, int RowNum) {
    for (int i = rank * RowNum; i < (rank + 1) * RowNum; i++) {
        X[i] = Xnew[i];
    }
}

void solver::devide(std::vector<double>& X, std::vector<std::vector<double>> A, int rank, int RowNum) {
    for (int i = rank * RowNum; i < (rank + 1) * RowNum; i++) {
        X[i] /= A[i][i];
    }
}

void solver::assignment(std::vector<double>& X, std::vector<double> F, int rank, int RowNum) {
    for (int i = rank * RowNum; i < (rank + 1) * RowNum; i++) {
        X[i] = F[i];
    }
}
