#include <iostream>
#include <chrono>
#include <random>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <map>

// ---------------- Matrix helpers ----------------
double** createMatrix(int n) {
    double** matrix = new double*[n];
    for (int i = 0; i < n; ++i) matrix[i] = new double[n];
    return matrix;
}

void deleteMatrix(double** matrix, int n) {
    for (int i = 0; i < n; ++i) delete[] matrix[i];
    delete[] matrix;
}

void initializeMatrices(double** A, double** B, double** C, int n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            B[i][j] = dis(gen);
            C[i][j] = dis(gen);
            A[i][j] = 0.0;
        }
    }
}

// ---------------- Multiplication orders ----------------
void multiply_ijk(double** A, double** B, double** C, int n) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            A[i][j] = 0.0;
            for (int k = 0; k < n; ++k)
                A[i][j] += B[i][k] * C[k][j];
        }
}

void multiply_ikj(double** A, double** B, double** C, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) A[i][j] = 0.0;
        for (int k = 0; k < n; ++k)
            for (int j = 0; j < n; ++j)
                A[i][j] += B[i][k] * C[k][j];
    }
}

void multiply_kij(double** A, double** B, double** C, int n) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) A[i][j] = 0.0;
    for (int k = 0; k < n; ++k)
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                A[i][j] += B[i][k] * C[k][j];
}

void multiply_kji(double** A, double** B, double** C, int n) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) A[i][j] = 0.0;
    for (int k = 0; k < n; ++k)
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i)
                A[i][j] += B[i][k] * C[k][j];
}

void multiply_jik(double** A, double** B, double** C, int n) {
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            A[i][j] = 0.0;
            for (int k = 0; k < n; ++k)
                A[i][j] += B[i][k] * C[k][j];
        }
}

void multiply_jki(double** A, double** B, double** C, int n) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) A[i][j] = 0.0;
    for (int j = 0; j < n; ++j)
        for (int k = 0; k < n; ++k)
            for (int i = 0; i < n; ++i)
                A[i][j] += B[i][k] * C[k][j];
}

// ---------------- Timing ----------------
double measureTime(void (*multiplyFunc)(double**, double**, double**, int),
                   double** A, double** B, double** C, int n) {
    auto start = std::chrono::high_resolution_clock::now();
    multiplyFunc(A, B, C, n);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    return duration.count(); // seconds
}

// Format helper: convert dot to comma for RU Excel
std::string fmt(double x, int prec = 6) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(prec) << x;
    std::string s = ss.str();
    std::replace(s.begin(), s.end(), '.', ',');
    return s;
}

// --------- Least squares cubic fit: y ≈ a3*x^3 + a2*x^2 + a1*x + a0 ---------
struct Poly3 { double a3, a2, a1, a0; };

static bool solve4(double M[4][4], double b[4], double x[4]) {
    int n = 4;
    for (int i=0;i<n;i++) {
        int br = i;
        double mv = std::fabs(M[i][i]);
        for (int r=i+1;r<n;r++) {
            double v = std::fabs(M[r][i]);
            if (v > mv) { mv = v; br = r; }
        }
        if (mv < 1e-18) return false;
        if (br != i) {
            for (int c=i;c<n;c++) std::swap(M[i][c], M[br][c]);
            std::swap(b[i], b[br]);
        }
        double diag = M[i][i];
        for (int c=i;c<n;c++) M[i][c] /= diag;
        b[i] /= diag;
        for (int r=0;r<n;r++) if (r!=i) {
            double f = M[r][i];
            if (std::fabs(f) < 1e-18) continue;
            for (int c=i;c<n;c++) M[r][c] -= f * M[i][c];
            b[r] -= f * b[i];
        }
    }
    for (int i=0;i<n;i++) x[i] = b[i];
    return true;
}

Poly3 polyfit3(const std::vector<double>& x, const std::vector<double>& y) {
    double S[7] = {0};
    for (double xv : x) {
        double p = 1.0;
        S[0] += 1.0;
        for (int k=1;k<=6;k++) { p *= xv; S[k] += p; }
    }
    double M[4][4];
    for (int p=0;p<4;p++)
        for (int q=0;q<4;q++)
            M[p][q] = S[p+q];
    double b[4] = {0,0,0,0};
    for (size_t i=0;i<x.size();++i) {
        double xp[4] = {1.0, x[i], x[i]*x[i], x[i]*x[i]*x[i]};
        for (int p=0;p<4;p++) b[p] += y[i]*xp[p];
    }
    double a[4];
    if (!solve4(M,b,a)) return {0,0,0,0};
    return {a[3], a[2], a[1], a[0]};
}

double r2_poly3(const std::vector<double>& x, const std::vector<double>& y, const Poly3& P) {
    double ymean = 0.0;
    for (double v: y) ymean += v;
    ymean /= (double)y.size();
    auto f = [&](double t){ return ((P.a3*t + P.a2)*t + P.a1)*t + P.a0; };
    double ss_tot=0.0, ss_res=0.0;
    for (size_t i=0;i<x.size();++i) {
        double yi = y[i];
        double fi = f(x[i]);
        ss_tot += (yi - ymean)*(yi - ymean);
        ss_res += (yi - fi)*(yi - fi);
    }
    if (ss_tot <= 1e-18) return 1.0;
    return 1.0 - ss_res/ss_tot;
}

void print_poly3(const std::string& name, const Poly3& P, double r2) {
    std::cout.setf(std::ios::scientific);
    std::cout << std::setprecision(6);
    std::cout << name << " trend: y = "
              << P.a3 << "*n^3 + "
              << P.a2 << "*n^2 + "
              << P.a1 << "*n + "
              << P.a0
              << " ; R^2 = " << r2 << "\n";
    std::cout.unsetf(std::ios::scientific);
}

// ---------------- Benchmark one size ----------------
struct Timings {
    double ijk, ikj, kij, kji, jik, jki;
};

Timings benchmarkForSize(int n, int repetitions, std::ofstream& csv) {
    std::cout << "\nMatrix size: " << n << "x" << n << "\n"
              << "==========================================\n";

    double** A = createMatrix(n);
    double** B = createMatrix(n);
    double** C = createMatrix(n);
    initializeMatrices(A, B, C, n);

    auto avg = [&](auto f){ double tot=0; for(int r=0;r<repetitions;++r) tot+=measureTime(f,A,B,C,n); return tot/repetitions; };

    const double avg_ijk = avg(multiply_ijk);
    const double avg_ikj = avg(multiply_ikj);
    const double avg_kij = avg(multiply_kij);
    const double avg_kji = avg(multiply_kji);
    const double avg_jik = avg(multiply_jik);
    const double avg_jki = avg(multiply_jki);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "(i,j,k): " << avg_ijk << "s\n";
    std::cout << "(i,k,j): " << avg_ikj << "s\n";
    std::cout << "(k,i,j): " << avg_kij << "s\n";
    std::cout << "(k,j,i): " << avg_kji << "s\n";
    std::cout << "(j,i,k): " << avg_jik << "s\n";
    std::cout << "(j,k,i): " << avg_jki << "s\n";

    auto gflops = [&](double t){ return (2.0 * n * n * n) / t / 1e9; };

    csv << "i,j,k;" << n << ";" << fmt(avg_ijk) << ";" << fmt(gflops(avg_ijk)) << "\n";
    csv << "i,k,j;" << n << ";" << fmt(avg_ikj) << ";" << fmt(gflops(avg_ikj)) << "\n";
    csv << "k,i,j;" << n << ";" << fmt(avg_kij) << ";" << fmt(gflops(avg_kij)) << "\n";
    csv << "k,j,i;" << n << ";" << fmt(avg_kji) << ";" << fmt(gflops(avg_kji)) << "\n";
    csv << "j,i,k;" << n << ";" << fmt(avg_jik) << ";" << fmt(gflops(avg_jik)) << "\n";
    csv << "j,k,i;" << n << ";" << fmt(avg_jki) << ";" << fmt(gflops(avg_jki)) << "\n";

    double min_time = avg_ijk; std::string best = "(i,j,k)";
    auto upd = [&](double t,const char* nstr){ if (t<min_time){ min_time=t; best=nstr; } };
    upd(avg_ikj,"(i,k,j)"); upd(avg_kij,"(k,i,j)"); upd(avg_kji,"(k,j,i)"); upd(avg_jik,"(j,i,k)"); upd(avg_jki,"(j,k,i)");
    std::cout << "Best method: " << best << " (" << min_time << "s)\n";

    deleteMatrix(A, n); deleteMatrix(B, n); deleteMatrix(C, n);
    return {avg_ijk, avg_ikj, avg_kij, avg_kji, avg_jik, avg_jki};
}

int main() {
    std::cout << "Matrix Multiplication Performance Comparison\n"
              << "==========================================\n";

    std::ofstream csv("results.csv");
    csv << "order;n;avg_time_sec;gflops\n";

    // ⬇⬇⬇ ДОБАВЛЕН 2500 и синхронизированы повторы
    const int sizes[]       = {100, 400, 700, 1000, 1500};
    const int repetitions[] = {5,   3,   2,    1,    1   };
    const int NUM_SIZES = static_cast<int>(sizeof(sizes)/sizeof(sizes[0]));

    std::vector<double> ns;
    ns.reserve(NUM_SIZES);

    std::map<std::string, std::vector<double>> ys;
    ys["i,j,k"] = {}; ys["i,k,j"] = {}; ys["k,i,j"] = {};
    ys["k,j,i"] = {}; ys["j,i,k"] = {}; ys["j,k,i"] = {};
    for (auto& kv : ys) kv.second.reserve(NUM_SIZES);

    for (int idx = 0; idx < NUM_SIZES; ++idx) {
        int n = sizes[idx];
        Timings t = benchmarkForSize(n, repetitions[idx], csv);
        ns.push_back(static_cast<double>(n));

        ys["i,j,k"].push_back(t.ijk);
        ys["i,k,j"].push_back(t.ikj);
        ys["k,i,j"].push_back(t.kij);
        ys["k,j,i"].push_back(t.kji);
        ys["j,i,k"].push_back(t.jik);
        ys["j,k,i"].push_back(t.jki);
    }

    csv.close();
    std::cout << "\nРезультаты сохранены в файл: results.csv\n";

    // ---------- Trend for each order ----------
    std::cout << "\n=== Polynomial trend (3rd degree) for each order ===\n";
    for (const auto& kv : ys) {
        const std::string& name = kv.first;
        const auto& y = kv.second;
        Poly3 P = polyfit3(ns, y);
        double r2 = r2_poly3(ns, y, P);
        print_poly3(name, P, r2);
    }

    // ---------- Overall average trend ----------
    std::vector<double> yavg(ns.size(), 0.0);
    for (size_t i = 0; i < ns.size(); ++i) {
        for (const auto& kv : ys) yavg[i] += kv.second[i];
        yavg[i] /= static_cast<double>(ys.size());
    }
    Poly3 Pavg = polyfit3(ns, yavg);
    double r2avg = r2_poly3(ns, yavg, Pavg);
    std::cout << "\nОбщий тренд (среднее по всем порядкам):\n";
    print_poly3("avg(all orders)", Pavg, r2avg);

    std::cout << "\n=== Analysis (short) ===\n"
              << "- Перестановка циклов не меняет математику, но меняет паттерн доступа к памяти.\n"
              << "- Порядки, которые проходят по строкам/столбцам последовательно, лучше используют кэш.\n";
    return 0;
}
