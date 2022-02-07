#include "Header.h"
#include <iostream>

double solver::df_bulk(int compIdx, int spatialIdx) {
    
    double n_sum = 0;
    double sum1 = 0, sum2 = 0;

    int M = settings["M"][0];
    double Temp = settings["Temp"][0];
    double b = settings["b"][0];

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
    int M = settings["M"][0]; 
    int N = settings["N"][0];
    for (size_t k = 0; k < M; k++) { //номер блока по горизонтали 
        for (size_t p = 0; p < M; p++) { // номер блока по вертикали 
            for (size_t i = k * (N); i < (1 + k) * N; i++) { // номер строки в этом блоке
                for (size_t j = p * (N); j < (1 + p) * N; j++) { // номер столбца в этом блоке

                    if ( (i == k*N) || (i == (k+1)*N - 1) ) { 
                        //A[i][j] = 0.0;
                        continue;
                    }
                    if ( (i - k * N) == (j - p * N) ) {
                        //A[i][j] = 1.0;
                        A.insert(std::make_pair(std::make_pair(i, j), 1.0));
                    }
                }
            }
        }
    }

    for (size_t k = 0; k < M; k++) { //номер блока по горизонтали 
        for (size_t p = 0; p < M; p++) { // номер блока по вертикали
            for (size_t i = k * N; i < (1 + k) * N; i++) { // номер строки в этом блоке
                for (size_t j = N * M + p * N; j < N * M + (1 + p) * N; j++) { // номер столбца в этом блоке
                    double coeff_1 = Mcoeff[k][p] * settings["delta_t"][0] / (settings["delta_z"][0] * settings["delta_z"][0]);
                    if ( (i == k * N) || (i == (k + 1) * N - 1) ) { 
                        //A[i][j] = 0.0;
                        continue;
                    }
                    if ( (i - k * N) == (j - N * M - p * N - 1) ) {
                        //A[i][j] = coeff_1;
                        A.insert(std::make_pair(std::make_pair(i, j), coeff_1));
                    }
                    if ( (i - k * N) == (j - N * M - p * N) ) {
                        //A[i][j] = -2*coeff_1;
                        A.insert(std::make_pair(std::make_pair(i, j), -2*coeff_1));
                    }
                    if ((i - k * N) == (j - N * M - p * N + 1)) {
                        //A[i][j] = coeff_1;
                        A.insert(std::make_pair(std::make_pair(i, j), coeff_1));
                    }
                }
            }
        }
    }
    
    for (size_t k = 0; k < M; k++) { //номер блока по горизонтали 
        for (size_t p = 0; p < M; p++) { // номер блока по вертикали 
            for (size_t i = N * M + k * N; i < N * M + (1 + k) * N; i++) { // номер строки в этом блоке
                for (size_t j = p * (N); j < (1 + p) * N; j++) { // номер столбца в этом блоке
                    double coeff_1 = c[k][p] / (settings["delta_z"][0] * settings["delta_z"][0]);
                    if ((i == N * M + k * N) || (i == (k + 1) * N - 1 + N * M)) {
                        //A[i][j] = 0.0;
                        continue;
                    }
                    if ((i - k * N - N * M) == (j  - p * N - 1)) {
                        //A[i][j] = coeff_1;
                        A.insert(std::make_pair(std::make_pair(i, j), coeff_1));
                    }
                    if ((i - k * N - N * M) == (j - p * N)) {
                        //A[i][j] = -2 * coeff_1;
                        A.insert(std::make_pair(std::make_pair(i, j), -2*coeff_1));
                    }
                    if ((i - k * N - N * M) == (j - p * N + 1)) {
                        //A[i][j] = coeff_1;
                        A.insert(std::make_pair(std::make_pair(i, j), coeff_1));
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
                        //A[i][j] = 0.0;
                        continue;
                    }
                    
                    if ((i - k * N - N * M) == (j - p * N - N * M)) {
                        //A[i][j] = 1.0;
                        A.insert(std::make_pair(std::make_pair(i, j), 1.0));
                    }else{ /*A[i][j] = 0.0;*/ }
                    
                }
            }
        }
    }
}

double solver::normaInf(std::vector<double> X, int rank) {  
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

void solver::mult_matvec(std::map<std::pair<int, int>, double> A, std::vector<double> b, std::vector<double>& Ab, int rank) {

    for (int i = rank * RowNum; i < (rank + 1) * RowNum; i++) {
        for (int j = 0; j < 2 * settings["N"][0] * settings["M"][0]; j++) {
            if (i != j)
                Ab[i] += A.find(std::make_pair(i, j))->second * b[j];
        }
    }
}

void solver::mult_matvec_(std::map<std::pair<int, int>, double> A, std::vector<double> b, std::vector<double>& Ab, int rank) {

    for (int i = rank * RowNum; i < (rank + 1) * RowNum; i++) {
        for (int j = 0; j < 2 * settings["N"][0] * settings["M"][0]; j++) {
            Ab[i] += A.find(std::make_pair(i, j))->second * b[j];
        }
    }
}

void solver::Jacobi(int rank) {
	
    std::vector<double> TempX;
    std::vector<double> Ax;
    std::vector<double> rezX;
	
	double norm = 0.0;
	int iter = 0;

	for (int i = 0; i < 2*settings["N"][0]*settings["M"][0]; i++) {
		Ax[i] = 0.0;
		TempX[i] = 0.0;
		rezX[i] = 0.0;
	}

	do {
		iter++;
		assignment(TempX, f, rank);
		mult_matvec(A, x, Ax, rank);
		subtract_vec(TempX, Ax, TempX, rank);
		devide(TempX, A, rank);

		abs_subtract_vec(x, TempX, rezX, rank);
		norm = normaInf(rezX, rank);
		if (rank == 0) { std::cout << "NORMA = " << norm << " " << "ITER " << iter << std::endl; }

		update_X(TempX, x, rank);

	} while (norm > settings["eps"][0]);
	if (rank == 0) { std::cout << "	Number of iteration = " << iter << std::endl; }
}

void solver::subtract_vec(std::vector<double> x, std::vector<double> y, std::vector<double>& rez, int rank) {
    for (int i = rank * RowNum; i < (rank + 1) * RowNum; i++) {
        rez[i] = x[i] - y[i];
    }
}

void solver::abs_subtract_vec(std::vector<double> x, std::vector<double> y, std::vector<double>& z, int rank) {
    for (int i = rank * RowNum; i < (rank + 1) * RowNum; i++) {
        z[i] = fabs(x[i] - y[i]);
    }
}

void solver::update_X(std::vector<double>& Xnew, std::vector<double>& X, int rank) {
    for (int i = rank * RowNum; i < (rank + 1) * RowNum; i++) {
        X[i] = Xnew[i];
    }
}

void solver::devide(std::vector<double>& X, std::map<std::pair<int, int>, double> A, int rank) {
    for (int i = rank * RowNum; i < (rank + 1) * RowNum; i++) {
        //X[i] /= A[i][i];
        X[i] /= A.find(std::make_pair(i, i))->second;
    }
}

void solver::assignment(std::vector<double>& X, std::vector<double> F, int rank) {
    for (int i = rank * RowNum; i < (rank + 1) * RowNum; i++) {
        X[i] = F[i];
    }
}

void solver::initialization(int world_size) {
    numOfProc = world_size;
    RowNum = 2*settings["N"][0]*settings["M"][0] / numOfProc;
    
    double omega = settings["omega"][0];
    if (omega <= 0.491) {
        m_omega = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
    }
    else {
        m_omega = 0.379642 + 1.48503 * omega - 0.164423 * omega * omega + 0.016666 * omega * omega * omega;
    }

    std::ifstream ifs("f:\\C++\\DensityGradientMethod\\Settings.txt");
    if (ifs.is_open()) {
        while (ifs) {
            std::string key, str;
            double v;
            std::vector<double> val;
            std::getline(ifs, key, ' ');
            std::getline(ifs, str, '\n');

            std::string delimiter = " ";

            size_t pos = 0;
            std::string token;
            while ((pos = str.find(delimiter)) != std::string::npos) {
                token = str.substr(0, pos);
                std::cout << token << std::endl;
                str.erase(0, pos + delimiter.length());
                std::stringstream s;
                s << token;
                s >> v;
                val.push_back(v);
            }
            std::stringstream s;
            s << str;
            s >> v;
            val.push_back(v);
            settings.emplace(key, val);
        }
    }
    
    for (int i = 0; i < settings["M"][0]; i++) {
        Tr.push_back(settings["Temp"][0] / settings["Tc"][i]);
        alpha_Tr_omega.push_back((1 + m_omega * (1 - sqrt(Tr[i]))) * (1 + m_omega * (1 - sqrt(Tr[i]))));
        b_i.push_back(0.0778 * settings["R"][0] * settings["Tc"][i] / settings["Pc"][i]);
        a_i.push_back(0.457235 * alpha_Tr_omega[i] * settings["R"][0] * settings["R"][0] * settings["Tc"][i] * settings["Tc"][i] / settings["Pc"][i]);
    }
    
    // input c, Mcoeff, Acoeff !
    
    RightSideVector();
    matrixA();
}
