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
    
    for (size_t i = 0; i < settings["M"][0]; i++) {
        for (size_t j = 0; j < settings["N"][0]; j++) {
            f.push_back(n[i][j]);
        }
    }

    for (size_t i = 0; i < settings["M"][0]; i++) {
        for (size_t j = 0; j < settings["N"][0]; j++) {
            f.push_back(df_bulk(i,j));
        }
    }
}

void solver::printA() {
    std::cout << "size = " << Arecv_value.size() << std::endl;
    for (int i = 0; i < Arecv_value.size(); i++) {
        std::cout << "A[" << Arecv_strIdx[i] << ", " << Arecv_colIdx[i] << "] = " << Arecv_value[i] << std::endl;
    }
}

void solver::matrixA() {
    int M = settings["M"][0]; 
    int N = settings["N"][0];
    for (size_t k = 0; k < M; k++) { //номер блока по горизонтали 
        for (size_t p = 0; p < M; p++) { // номер блока по вертикали 
            for (size_t i = k * N; i < (1 + k) * N; i++) { // номер строки в этом блоке
                for (size_t j = p * N; j < (1 + p) * N; j++) { // номер столбца в этом блоке
                    if (k == p) {
                        if ((i == k * N)) { // boundary condition on right side
                            if (i == j) {
                                A_strIdx.push_back(i);
                                A_colIdx.push_back(j);
                                A_value.push_back(1.0);
                                continue;
                            }
                        }

                        if ((i == (k + 1) * N - 1)) { // boundary condition on left side
                            if (i == j) {
                                A_strIdx.push_back(i);
                                A_colIdx.push_back(j);
                                A_value.push_back(1.0);
                                continue;
                            }
                        }
                        if ((i - k * N) == (j - p * N)) {
                            A_strIdx.push_back(i);
                            A_colIdx.push_back(j);
                            A_value.push_back(1.0);
                        }
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

                    if (coeff_1 == 0) { continue; }

                    if ((i == k * N)) { // boundary condition on right side
                        if (i - k * N == j - N * M - p * N) {
                            A_strIdx.push_back(i);
                            A_colIdx.push_back(j);
                            A_value.push_back(-coeff_1);
                            continue;
                        }
                        if (i - k * N == j - 1 - N * M - p * N) {
                            A_strIdx.push_back(i);
                            A_colIdx.push_back(j);
                            A_value.push_back(coeff_1);
                            continue;
                        }
                    }

                    if ((i == (k + 1) * N - 1)) { // boundary condition on left side
                        if (i - k * N == j - N * M - p * N) {
                            A_strIdx.push_back(i);
                            A_colIdx.push_back(j);
                            A_value.push_back(-coeff_1);
                            continue;
                        }
                        if (i - k * N == j - N * M - p * N - 1) {
                            A_strIdx.push_back(i);
                            A_colIdx.push_back(j);
                            A_value.push_back(coeff_1);
                            continue;
                        }
                    }

                    if ((i - k * N) == (j - N * M - p * N - 1)) {
                        A_strIdx.push_back(i);
                        A_colIdx.push_back(j);
                        A_value.push_back(coeff_1);
                    }
                    if ((i - k * N) == (j - N * M - p * N)) {
                        A_strIdx.push_back(i);
                        A_colIdx.push_back(j);
                        A_value.push_back(-2 * coeff_1);
                    }
                    if ((i - k * N) == (j - N * M - p * N + 1)) {
                        A_strIdx.push_back(i);
                        A_colIdx.push_back(j);
                        A_value.push_back(coeff_1);
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

                    if (coeff_1 == 0) { continue; }

                    if ((i == N * M + k * N)) { // boundary condition on right side
                        if (i - k * N - N * M == j - p * N) {
                            A_strIdx.push_back(i);
                            A_colIdx.push_back(j);
                            A_value.push_back(-coeff_1);
                            continue;
                        }
                        if (i - k * N - N * M == j - 1 - p * N) {
                            A_strIdx.push_back(i);
                            A_colIdx.push_back(j);
                            A_value.push_back(coeff_1);
                            continue;
                        }
                    }

                    if ((i == (k + 1) * N - 1 + N * M)) { // boundary condition on left side
                        if (i - k * N - N * M == j - p * N) {
                            A_strIdx.push_back(i);
                            A_colIdx.push_back(j);
                            A_value.push_back(-coeff_1);
                            continue;
                        }
                        if (i - k * N - N * M == j - p * N - 1) {
                            A_strIdx.push_back(i);
                            A_colIdx.push_back(j);
                            A_value.push_back(coeff_1);
                            continue;
                        }
                    }

                    if ((i - k * N - N * M) == (j - p * N - 1)) {
                        A_strIdx.push_back(i);
                        A_colIdx.push_back(j);
                        A_value.push_back(coeff_1);
                    }
                    if ((i - k * N - N * M) == (j - p * N)) {
                        A_strIdx.push_back(i);
                        A_colIdx.push_back(j);
                        A_value.push_back(-2 * coeff_1);
                    }
                    if ((i - k * N - N * M) == (j - p * N + 1)) {
                        A_strIdx.push_back(i);
                        A_colIdx.push_back(j);
                        A_value.push_back(coeff_1);
                    }
                }
            }
        }
    }

    for (size_t k = 0; k < M; k++) { //номер блока по горизонтали 
        for (size_t p = 0; p < M; p++) { // номер блока по вертикали 
            for (size_t i = N * M + k * N; i < N * M + (1 + k) * N; i++) { // номер строки в этом блоке
                for (size_t j = N * M + p * N; j < N * M + (1 + p) * N; j++) { // номер столбца в этом блоке
                    if (k == p) {
                        if ((i == N * M + k * N)) { // boundary condition on right side
                            if (i == j) {
                                A_strIdx.push_back(i);
                                A_colIdx.push_back(j);
                                A_value.push_back(1.0);
                                continue;
                            }
                        }

                        if ((i == (k + 1) * N - 1 + N * M)) { // boundary condition on left side
                            if (i == j) {
                                A_strIdx.push_back(i);
                                A_colIdx.push_back(j);
                                A_value.push_back(1.0);
                                continue;
                            }
                        }

                        if ((i - k * N - N * M) == (j - p * N - N * M)) {
                            A_strIdx.push_back(i);
                            A_colIdx.push_back(j);
                            A_value.push_back(1.0);
                        }
                    }
                }
            }
        }
    }

    sendcount = A_value.size() / numOfProc;

    Arecv_value.resize(sendcount);
    Arecv_colIdx.resize(sendcount);
    Arecv_strIdx.resize(sendcount);

    MPI_Scatter(A_value.data(), sendcount, MPI_DOUBLE, Arecv_value.data(), sendcount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(A_strIdx.data(), sendcount, MPI_INT, Arecv_colIdx.data(), sendcount, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(A_colIdx.data(), sendcount, MPI_INT, Arecv_strIdx.data(), sendcount, MPI_INT, 0, MPI_COMM_WORLD);
}

double solver::normaInf() {  
    double norm = 0.0;   
    norm = fabs(TempX[0]);
    for (int i = 0; i < TempX.size(); i++) {
        if (fabs(TempX[i]) > norm) {
            norm = fabs(TempX[i]);
        }
    }
    return norm;
}

void solver::mult_matvec(int rank) {
    for (int i = 0; i < sendcount; i++) {
        if (Arecv_colIdx[i] != Arecv_strIdx[i]) {
            Ax[i+rank*sendcount] += Arecv_value[i] * x[i + rank * sendcount];
        }
    }
}

void solver::Jacobi(int rank) {
    Ax.resize(2 * settings["N"][0] * settings["M"][0]);
	double norm = 0.0;
	int iter = 0;

	do {
		iter++;
		assignment();
		mult_matvec(rank);
		subtract_vec();
		devide();

		abs_subtract_vec();
		norm = normaInf();
		if (rank == 0) { std::cout << "NORMA = " << norm << " " << "ITER " << iter << std::endl; }

		update_X();

	} while (norm > settings["eps"][0]);
	if (rank == 0) { std::cout << "	Number of iteration = " << iter << std::endl; }
}

void solver::assignment() {
    for (int i = 0; i < f.size(); i++) { //int i = rank * sendcount; i < (rank + 1) * sendcount; i++
        TempX.push_back(f[i]);
    }
}

void solver::subtract_vec() {
    for (int i = 0; i < TempX.size(); i++) {
        TempX[i] = TempX[i] - Ax[i];
    }
}

void solver::abs_subtract_vec() {
    for (int i = 0; i < TempX.size(); i++) {
        TempX[i] = fabs(x[i] - TempX[i]);
    }
}

void solver::update_X() {
    x.resize(TempX.size());
    for (int i = 0; i < TempX.size(); i++) {
        x[i] = TempX[i];
    }
}

void solver::devide() {
    
    for (int i = 0; i < sendcount; i++) {
        if (Arecv_strIdx[i] == Arecv_colIdx[i]) {
            TempX[Arecv_strIdx[i]] /=  Arecv_value[i];
        }
    }
}

void solver::initialization(int world_size) {
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
        b_i.push_back(0.0778 * R * settings["Tc"][i] / settings["Pc"][i]);
        a_i.push_back(0.457235 * alpha_Tr_omega[i] * R * R * settings["Tc"][i] * settings["Tc"][i] / settings["Pc"][i]);
    }

    numOfProc = world_size;

    double omega = settings["omega"][0];
    if (omega <= 0.491) {
        m_omega = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
    }
    else {
        m_omega = 0.379642 + 1.48503 * omega - 0.164423 * omega * omega + 0.016666 * omega * omega * omega;
    }

    Mcoeff_init();

    Acoeff_init();
   
    Ccoeff_init();
    
    n_initDist_init();

    x.resize(2*settings["M"][0] * settings["N"][0]);

    RightSideVector();

    matrixA();

}

void solver::Mcoeff_init() {
    std::ifstream ifs("f:\\C++\\DensityGradientMethod\\Mcoeff.txt");
    Mcoeff.resize(settings["M"][0]);
    if (ifs.is_open()) {
        int i = 0;
        while (ifs) {
            std::string str;
            double m;
            std::getline(ifs, str, '\n');
            std::string delimiter = " ";

            size_t pos = 0;
            std::string token;
            while ((pos = str.find(delimiter)) != std::string::npos) {
                token = str.substr(0, pos);
                str.erase(0, pos + delimiter.length());
                std::stringstream s;
                s << token;
                s >> m;
                Mcoeff[i].push_back(m);
            }
            if (str == "") { break; }
            std::stringstream s;
            s << str;
            s >> m;
            Mcoeff[i].push_back(m);
            i++;
        }
    }
}

void solver::Acoeff_init() {
    std::ifstream ifs("f:\\C++\\DensityGradientMethod\\Acoeff.txt");
    Acoeff.resize(settings["M"][0]);
    if (ifs.is_open()) {
        int i = 0;
        while (ifs) {
            std::string str;
            double m;
            std::getline(ifs, str, '\n');
            std::string delimiter = " ";

            size_t pos = 0;
            std::string token;
            while ((pos = str.find(delimiter)) != std::string::npos) {
                token = str.substr(0, pos);
                str.erase(0, pos + delimiter.length());
                std::stringstream s;
                s << token;
                s >> m;
                Acoeff[i].push_back(m);
            }
            if (str == "") { break; }
            std::stringstream s;
            s << str;
            s >> m;
            Acoeff[i].push_back(m);
            i++;
        }
    }
}

void solver::Ccoeff_init() {
    std::ifstream ifs("f:\\C++\\DensityGradientMethod\\Ccoeff.txt");
    c.resize(settings["M"][0]);
    if (ifs.is_open()) {
        int i = 0;
        while (ifs) {
            std::string str;
            double m;
            std::getline(ifs, str, '\n');
            std::string delimiter = " ";

            size_t pos = 0;
            std::string token;
            while ((pos = str.find(delimiter)) != std::string::npos) {
                token = str.substr(0, pos);
                str.erase(0, pos + delimiter.length());
                std::stringstream s;
                s << token;
                s >> m;
                c[i].push_back(m);
            }
            if (str == "") { break; }
            std::stringstream s;
            s << str;
            s >> m;
            c[i].push_back(m);
            i++;
        }
    }
}

void solver::n_initDist_init() {
    std::ifstream ifs("f:\\C++\\DensityGradientMethod\\n_initial_distribution.txt");
    int M = settings["M"][0];
    n.resize(M);
    if (ifs.is_open()) {
        int i = 0;
        while (ifs) {
            std::string str;
            double m;
            std::getline(ifs, str, '\n');
            std::string delimiter = " ";

            size_t pos = 0;
            std::string token;
            while ((pos = str.find(delimiter)) != std::string::npos) {
                token = str.substr(0, pos);
                str.erase(0, pos + delimiter.length());
                std::stringstream s;
                s << token;
                s >> m;
                n[i].push_back(m);
            }
            if (str == "") { break; }
            std::stringstream s;
            s << str;
            s >> m;
            n[i].push_back(m);
            i++;
        }
    }
}