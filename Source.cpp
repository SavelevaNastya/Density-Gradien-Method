#include "Header.h"
#include <iostream>

double solver::df_bulk(int compIdx, int timeIdx, int spatialIdx) {
    double n_sum = 0;
    double sum1 = 0, sum2 = 0;
    for (int i = 0; i < M; i++) {
        n_sum += n[i][timeIdx][spatialIdx];
        sum1 += n[i][timeIdx][spatialIdx] * Acoeff[compIdx][i];
        sum2 += n[i][timeIdx][spatialIdx] * Acoeff[i][compIdx];
    }

    double rez = R * Temp * (b_i[compIdx] * n_sum / (1 - b * n_sum) + log(n[compIdx][timeIdx][spatialIdx]) - log(1 - b * n_sum))
        + (sum1 + sum2) * log((b * delta_1 * n_sum + 1) / (b * delta_2 * n_sum + 1)) / (delta_2 - delta_1)
        - a * b_i[compIdx] * n_sum / (b * (delta_1 * b * n_sum + 1)*(b * (delta_2 * b * n_sum + 1)));
    return rez;
}

double solver::f(int compIdx, int timeIdx, int spatialIdx) {
    double sum1 = 0;
    double sum2 = 0;
    for (int j = 0; j < M; j++) {
        for (int q = 0; q < M; q++) {
            sum2 += c[j][q] * (n[q][timeIdx][spatialIdx + 2] - 2 * n[q][timeIdx][spatialIdx + 1] + n[q][timeIdx][spatialIdx]
                + 2 * (n[q][timeIdx][spatialIdx + 1] - 2 * n[q][timeIdx][spatialIdx] + n[q][timeIdx][spatialIdx - 1])
                - n[q][timeIdx][spatialIdx] - 2 * n[q][timeIdx][spatialIdx - 1] + n[q][timeIdx][spatialIdx - 2])
                / (delta_z * delta_z);

        }
        sum1 += Mcoeff[compIdx][j] *
            (df_bulk(j, timeIdx - 1, spatialIdx + 1) + df_bulk(j, timeIdx - 1, spatialIdx - 1) - 2 * df_bulk(j, timeIdx - 1, spatialIdx)
                + sum2) / (delta_z * delta_z);
    }
    double rez = (n[compIdx][timeIdx][spatialIdx] - n[compIdx][timeIdx - 1][spatialIdx]) / delta_t + sum1;
    return rez;
}

void solver::init_cond() {
    define_m_omega();
}
void solver::FvectorFormation(int timeIdx, int spatialIdx) {
    for (int i = 0; i < M; i++) {
        F[i] = f(i, timeIdx, spatialIdx);
    }
}
void solver::WmatrixFormation(int compIdx, int timeIdx, int spatialIdx) { //Jacobi Matrix
   
}

void solver::matInverse(int compIdx, int timeIdx, int spatialIdx) { //inverse Jacobi matrix
    
};

double* solver::delta(int compIdx, int timeIdx, int spatialIdx) {
    //W_inv(n[compIdx, timeIdx, spatialIdx]) * F(n[compIdx, timeIdx, spatialIdx]);
    double rez[M];
    FvectorFormation(timeIdx, spatialIdx);
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            rez[i] = W_inv[i][j] * F[j];
        }
    }
    return rez;
}

double solver::norm(int timeIdx, int spatialIdx) { // || n_new - n || 
    double abs_vec[M];
    for (int compIdx = 0; compIdx < M; compIdx++) {
        abs_vec[compIdx] = n[compIdx][timeIdx][spatialIdx] - n_new[compIdx][timeIdx][spatialIdx];
     }
    double rez = fabs(abs_vec[0]);
    for (int i = 1; i < M; i++) {
        if (fabs(abs_vec[i]) > rez) {
            rez = fabs(abs_vec[i]);
        }
    }
    return rez;
}

void solver::newtonMethod(int timeIdx) {
    
    double norma = 0;
    int iter = 0;
    double* deltaPtr;
    
    for (int spatialIdx = 0; spatialIdx < M; spatialIdx++) { // последовательность циклов???
        while (norma > eps) {
            iter++;
            for (int compIdx = 0; compIdx < M; compIdx++) {
            deltaPtr = delta(compIdx, timeIdx, spatialIdx);
            n_new[compIdx][timeIdx][spatialIdx] = n[compIdx][timeIdx][spatialIdx] - *(deltaPtr + compIdx);
            norma = norm(timeIdx, spatialIdx);
            }
        }
    }
}

void solver::solve() {
    for (int timeIdx = 0; timeIdx < T; timeIdx = timeIdx + delta_t){
        newtonMethod(timeIdx);
    }
}