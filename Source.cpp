#include "Header.h"

using namespace Eigen;
#include <math.h>

double solver::df_bulk(int compIdx, int spatialIdx) {
    
    double n_sum = 0;
    double sum1 = 0, sum2 = 0;

    double Temp = settings["Temp"][0];
    double b = 0, a = 0;

    Acoeff.resize(M);
    for (int i = 0; i < M; i++) {
        Acoeff[i].resize(M);
    }
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            Acoeff[i][j] = sqrt(a_i[i] * a_i[j]) * (1 - Kcoeff[i][j]);
        }
    }
    std::vector<double> L_pure_comp;
    L_pure_comp.push_back(2.6147e-20);
    L_pure_comp.push_back(1.2415e-18);
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            c[i][j] = sqrt(L_pure_comp[i] * L_pure_comp[j]) * (1 - c[i][j]);
        }
    }
    for (int i = 0; i < M; i++) {
        n_sum += n[i][spatialIdx];
        sum1 += n[i][spatialIdx] * Acoeff[compIdx][i];
        sum2 += n[i][spatialIdx] * Acoeff[i][compIdx];
    }

    z.resize(M);
    for (int i = 0; i < M; i++) {
        z[i] = n[i][spatialIdx] / n_sum;
        b += z[i] * b_i[i];
    }

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            a += z[i] * z[j] * Acoeff[i][j];
        }
    }

    double rez_1 = R * Temp * (b_i[compIdx] * n_sum / (1 - b * n_sum) + log(n[compIdx][spatialIdx]) - log(1 - b * n_sum));
    double rez_2 = ((sum1 + sum2) / (b * n_sum) - a * b_i[compIdx] / (b * b)) * log((b * delta_1 * n_sum + 1) / (b * delta_2 * n_sum + 1)) / (delta_2 - delta_1);
    double rez_3 = a * b_i[compIdx] * n_sum / (b * (delta_1 * b * n_sum + 1) * (delta_2 * b * n_sum + 1));
    double rez = rez_1 + rez_2 - rez_3;
    return rez;
}

void solver::RightSideVector() {
    
    f.resize(2 * N * M);
    int count = 0;
    for (size_t i = 0; i < M; i++) {
        for (size_t j = 0; j < N; j++) {
            f(count) = n[i][j]; 
            count++;
        }
    }

    for (size_t i = 0; i < M; i++) {
        for (size_t j = 0; j < N; j++) {
            f(count) = df_bulk(i, j);
            count++;
        }
    }
}

void solver::matrixA() {
    
    int sz = 2 * M * N;
    A_value.resize(sz, sz);
    std::cout << "Initialization A: start..." << std::endl;
    for (size_t k = 0; k < M; k++) { //номер блока по горизонтали 
        std::cout << "Initialization A: first " << k << "..." << std::endl;
        for (size_t p = 0; p < M; p++) { // номер блока по вертикали 
            for (size_t i = k * N; i < (1 + k) * N; i++) { // номер строки в этом блоке
                for (size_t j = p * N; j < (1 + p) * N; j++) { // номер столбца в этом блоке
                    if (k == p) {
                        if ((i == k * N)) { // boundary condition on right side
                            if (i == j) {
                                A_value.insert(i, j) = 1.0;
                                continue;
                            }
                        }

                        if ((i == (k + 1) * N - 1)) { // boundary condition on left side
                            if (i == j) {
                                A_value.insert(i, j) = 1.0;
                                continue;
                            }
                        }
                        if ((i - k * N) == (j - p * N)) {
                            A_value.insert(i, j) = 1.0;
                        }
                    }
                }
            }
        }
    }
    std::cout << "  Initialization A, second block: start..." << std::endl;
    for (size_t k = 0; k < M; k++) { //номер блока по горизонтали 
        std::cout << "Initialization A: second " << k << "..." << std::endl;
        for (size_t p = 0; p < M; p++) { // номер блока по вертикали
            for (size_t i = k * N; i < (1 + k) * N; i++) { // номер строки в этом блоке
                for (size_t j = N * M + p * N; j < N * M + (1 + p) * N; j++) { // номер столбца в этом блоке
                    double coeff_1 = Mcoeff[k][p] * delta_t / (delta_z * delta_z);

                    if (coeff_1 == 0) { continue; }

                    if ((i == k * N)) { // boundary condition on right side
                        if (i - k * N == j - N * M - p * N) {
                            A_value.insert(i, j) = -coeff_1;
                            continue;
                        }
                        if (i - k * N == j - 1 - N * M - p * N) {
                            A_value.insert(i, j) = coeff_1;
                            continue;
                        }
                    }

                    if ((i == (k + 1) * N - 1)) { // boundary condition on left side
                        if (i - k * N == j - N * M - p * N) {
                            A_value.insert(i, j) = -coeff_1;
                            continue;
                        }
                        if (i - k * N == j - N * M - p * N - 1) {
                            A_value.insert(i, j) = coeff_1;
                            continue;
                        }
                    }

                    if ((i - k * N) == (j - N * M - p * N - 1)) {
                        A_value.insert(i, j) = coeff_1;
                    }
                    if ((i - k * N) == (j - N * M - p * N)) {
                        A_value.insert(i, j) = -2 * coeff_1;
                    }
                    if ((i - k * N) == (j - N * M - p * N + 1)) {
                        A_value.insert(i, j) = coeff_1;
                    }
                }
            }
        }
    }
    std::cout << "  Initialization A, third block: start..." << std::endl;
    for (size_t k = 0; k < M; k++) { //номер блока по горизонтали 
        std::cout << "Initialization A: third " << k << "..." << std::endl;
        for (size_t p = 0; p < M; p++) { // номер блока по вертикали 
            for (size_t i = N * M + k * N; i < N * M + (1 + k) * N; i++) { // номер строки в этом блоке
                for (size_t j = p * (N); j < (1 + p) * N; j++) { // номер столбца в этом блоке
                    double coeff_1 = c[k][p] / (delta_z * delta_z);

                    if (coeff_1 == 0) { continue; }

                    if ((i == N * M + k * N)) { // boundary condition on right side
                        if (i - k * N - N * M == j - p * N) {
                            A_value.insert(i, j) = -coeff_1;
                            continue;
                        }
                        if (i - k * N - N * M == j - 1 - p * N) {
                            A_value.insert(i, j) = coeff_1;
                            continue;
                        }
                    }

                    if ((i == (k + 1) * N - 1 + N * M)) { // boundary condition on left side
                        if (i - k * N - N * M == j - p * N) {
                            A_value.insert(i, j) = -coeff_1;
                            continue;
                        }
                        if (i - k * N - N * M == j - p * N - 1) {
                            A_value.insert(i, j) = coeff_1;
                            continue;
                        }
                    }

                    if ((i - k * N - N * M) == (j - p * N - 1)) {
                        A_value.insert(i, j) = coeff_1;
                    }
                    if ((i - k * N - N * M) == (j - p * N)) {
                        A_value.insert(i, j) = -2 * coeff_1;
                        ;
                    }
                    if ((i - k * N - N * M) == (j - p * N + 1)) {
                        A_value.insert(i, j) = coeff_1;
                    }
                }
            }
        }
    }
    std::cout << "  Initialization A, fourth block: start..." << std::endl;
    for (size_t k = 0; k < M; k++) { //номер блока по горизонтали 
        std::cout << "Initialization A: fourth " << k << "..." << std::endl;
        for (size_t p = 0; p < M; p++) { // номер блока по вертикали 
            for (size_t i = N * M + k * N; i < N * M + (1 + k) * N; i++) { // номер строки в этом блоке
                for (size_t j = N * M + p * N; j < N * M + (1 + p) * N; j++) { // номер столбца в этом блоке
                    if (k == p) {
                        if ((i == N * M + k * N)) { // boundary condition on right side
                            if (i == j) {
                                A_value.insert(i, j) = 1.0;
                                continue;
                            }
                        }

                        if ((i == (k + 1) * N - 1 + N * M)) { // boundary condition on left side
                            if (i == j) {
                                A_value.insert(i, j) = 1.0;
                                continue;
                            }
                        }

                        if ((i - k * N - N * M) == (j - p * N - N * M)) {
                            A_value.insert(i, j) = 1.0;
                        }
                    }
                }
            }
        }
    }
    std::cout << "Initialization A: done" << std::endl;
}

void solver::initialization() {
    
    std::ifstream ifs("./Parameters/Settings.txt");
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

    N = settings["D"][0] / settings["delta_z"][0];
    M = settings["M"][0];
    x.resize(2 * M * N);

    delta_z = settings["delta_z"][0];
    delta_t = settings["delta_t"][0];
    m_omega.resize(M);
    for (int i = 0; i < M; i++) {
        double omega = settings["omega"][i];
        if (omega <= 0.491) {
                m_omega[i] = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
        }
        else {
                m_omega[i] = 0.379642 + 1.48503 * omega - 0.164423 * omega * omega + 0.016666 * omega * omega * omega;
        }
    }

    for (int i = 0; i < M; i++) {
        Tr.push_back(settings["Temp"][0] / settings["Tc"][i]);
        alpha_Tr_omega.push_back((1 + m_omega[i] * (1 - sqrt(Tr[i]))) * (1 + m_omega[i] * (1 - sqrt(Tr[i]))));
        b_i.push_back(0.077796 * R * settings["Tc"][i] / settings["Pc"][i]);
        a_i.push_back(0.457235 * alpha_Tr_omega[i] * R * R * settings["Tc"][i] * settings["Tc"][i] / settings["Pc"][i]);
    }
    
    Mcoeff_init();

    Kcoeff_init();
   
    Lcoeff_init();
    
    n_initDist_init();

    RightSideVector();

    matrixA();
}

void solver::fill_matrix_from_file(std::string path, std::vector<std::vector<double>>& Matrix) {
    std::ifstream ifs(path);
    Matrix.resize(M);
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
                Matrix[i].push_back(m);
            }
            if (str == "") { break; }
            std::stringstream s;
            s << str;
            s >> m;
            Matrix[i].push_back(m);
            i++;
        }
    }
}

void solver::Mcoeff_init() {

    fill_matrix_from_file("./Parameters/Mcoeff.txt", Mcoeff);
}

void solver::Kcoeff_init() {

    fill_matrix_from_file("./Parameters/Kcoeff.txt", Kcoeff);
    
}

void solver::Lcoeff_init() {

    fill_matrix_from_file("./Parameters/Lcoeff.txt", c);
   
}

void solver::n_initDist_init() {
    int i;

    n.resize(M);
    for (i = 0; i < M; i++) { n[i].resize(N); }

    for (i = 0; i < N; i++) {
        if (i < N / 2) {
            n[0][i] = 200.0;
        }
        if (i >= N / 2) {
            n[0][i] = 99.0;
        }
    }
    for (i = 0; i < N; i++) {
        if (i < N / 2) {
            n[1][i] = 557.0;
        }
        if (i >= N / 2) {
            n[1][i] = 0.32;
        }
    }    
}

void solver::writeAnswer() {
    _mkdir("Results");

    _mkdir("Results/n");
    _mkdir("Results/mu");

    _mkdir("Results/n/CO2");
    _mkdir("Results/n/n_decane");
    _mkdir("Results/mu/CO2");
    _mkdir("Results/mu/n_decane");

    std::ostringstream strs;
    strs << time;
    std::string str = strs.str();

    std::string name_n_co2 = "Results/n/CO2/" + str + ".txt";
    std::ofstream n_co2(name_n_co2);

    std::string name_n_ndecane = "Results/n/n_decane/" + str + ".txt";
    std::ofstream n_ndecane(name_n_ndecane);

    std::string name_mu_co2 = "Results/mu/CO2/" + str + ".txt";
    std::ofstream mu_co2(name_mu_co2);

    std::string name_mu_ndecane = "Results/mu/n_decane/" + str + ".txt";
    std::ofstream mu_ndecane(name_mu_ndecane);

    

    for (int i = 0; i < x.size(); i++) {

        if(i < N){ n_co2 << x[i] << std::endl; }

        if (i > N && i < N * M) { n_ndecane << x[i] << std::endl; }

        if (i >= N * M && i < (1 + M) * N) { mu_co2 << x[i] << std::endl; }

        if (i >= (1 + M) * N) { mu_ndecane << x[i] << std::endl; }
    }

}

void solver::solve() {
    for (int t = 0; t < 100; t++){

        SimplicialCholesky<SparseMatrix<double>> chol(A_value);
        x = chol.solve(f);

        time = delta_t * t;
        std::cout << "  TIME = " << time << std::endl;
        writeAnswer();
        sigma();
        new_time_step();
    }
};

void solver::new_time_step() {
    
    int count = 0;
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            n[i][j] = x[count];
            count++;
        }
    }

    for (int i = 0; i < N * M; i++) {
        f[i] = x[i];
    }
    count = N * M;
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            f[count] = df_bulk(i, j);
            count++;
        }
    }
}

void solver::sigma() {
    double sigma = 0, sum = 0;
    
    int count = 0;
    n_new.resize(M);
    for (int i = 0; i < M; i++) { n_new[i].resize(N); }
    for(int i = 0; i < M; i++){
        for (int j = 0; j < N; j++) {
            n_new[i][j] = x[count];
            count++;
        }
    }

    for (int i = 1; i < N - 2; i++) {
        sum += f_sigma(i - 1) + 4 * f_sigma(i) + f_sigma(i + 1);
    }
    sigma = delta_z * sum / 3.0;
    std::cout << "      sigma = " << sigma << std::endl;
}

double solver::f_sigma(int spatialidx) {
    double sum = 0;
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            sum += c[i][j] * ((n_new[i][spatialidx + 1] - n_new[i][spatialidx]) / delta_z) * ((n_new[j][spatialidx + 1] - n_new[j][spatialidx]) / delta_z);
        }
    }
    
    return sum;
}