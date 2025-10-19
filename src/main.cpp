// ===============================
// Семинар 5 — допзадание
// Сравнение порядков трёх вложенных циклов при умножении матриц
// Генерируем CSV с временами, печатаем лучший порядок,
// считаем полиномиальную аппроксимацию (deg=3) и R^2 для каждого порядка.
// ===============================

#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <algorithm>

// -------------------------------
// Удобное короткое имя типа для матрицы n×n,
// храним как один сплошной массив (улучшает локальность кэша)
// Доступ: A[i*n + j]
// -------------------------------
using Matrix = std::vector<double>;

// Возврат смещения (индекса) одного элемента в сплошном массиве
inline size_t idx(size_t n, size_t i, size_t j) { return i * n + j; }

// -------------------------------
// Генерация случайных матриц (равномерное распределение [0;1))
// Фиксируем seed для воспроизводимости
// -------------------------------
void fill_random(Matrix& M, size_t n, uint32_t seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (size_t i = 0; i < n * n; ++i) M[i] = dist(rng);
}

// -------------------------------
// Таймер на std::chrono для замеров
// -------------------------------
struct Timer {
    using clock = std::chrono::high_resolution_clock;
    clock::time_point t0;
    void start() { t0 = clock::now(); }
    double stop_sec() const {
        auto dt = clock::now() - t0;
        return std::chrono::duration<double>(dt).count();
    }
};

// -------------------------------
// Умножение A = B * C разными порядками циклов
// Все функции предполагают нулевую инициализацию A перед вызовом
// -------------------------------

// Порядок (i, j, k)
void mul_ijk(Matrix& A, const Matrix& B, const Matrix& C, size_t n) {
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j) {
            double sum = 0.0;                 // локальный аккумулятор быстрее, чем писать A[...] += ...
            for (size_t k = 0; k < n; ++k)
                sum += B[idx(n,i,k)] * C[idx(n,k,j)];
            A[idx(n,i,j)] = sum;
        }
}

// Порядок (i, k, j)
void mul_ikj(Matrix& A, const Matrix& B, const Matrix& C, size_t n) {
    for (size_t i = 0; i < n; ++i)
        for (size_t k = 0; k < n; ++k) {
            double bik = B[idx(n,i,k)];
            for (size_t j = 0; j < n; ++j)
                A[idx(n,i,j)] += bik * C[idx(n,k,j)]; // выгодно: j идёт подряд → хороший проход по строке A и C
        }
}

// Порядок (k, i, j)
void mul_kij(Matrix& A, const Matrix& B, const Matrix& C, size_t n) {
    for (size_t k = 0; k < n; ++k)
        for (size_t i = 0; i < n; ++i) {
            double bik = B[idx(n,i,k)];
            for (size_t j = 0; j < n; ++j)
                A[idx(n,i,j)] += bik * C[idx(n,k,j)];
        }
}

// Порядок (k, j, i)
void mul_kji(Matrix& A, const Matrix& B, const Matrix& C, size_t n) {
    for (size_t k = 0; k < n; ++k)
        for (size_t j = 0; j < n; ++j) {
            double ckj = C[idx(n,k,j)];
            for (size_t i = 0; i < n; ++i)
                A[idx(n,i,j)] += B[idx(n,i,k)] * ckj;
        }
}

// Порядок (j, i, k)
void mul_jik(Matrix& A, const Matrix& B, const Matrix& C, size_t n) {
    for (size_t j = 0; j < n; ++j)
        for (size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            for (size_t k = 0; k < n; ++k)
                sum += B[idx(n,i,k)] * C[idx(n,k,j)];
            A[idx(n,i,j)] = sum;
        }
}

// Порядок (j, k, i)
void mul_jki(Matrix& A, const Matrix& B, const Matrix& C, size_t n) {
    for (size_t j = 0; j < n; ++j)
        for (size_t k = 0; k < n; ++k) {
            double ckj = C[idx(n,k,j)];
            for (size_t i = 0; i < n; ++i)
                A[idx(n,i,j)] += B[idx(n,i,k)] * ckj;
        }
}

// -------------------------------
// Измерение времени одной функции умножения
// repeats — число повторов (усредняем для стабильности)
// -------------------------------
template <typename MulFn>
double measure_mul(MulFn mul, size_t n, unsigned repeats,
                   const Matrix& B, const Matrix& C, Matrix& A)
{
    double acc = 0.0;
    for (unsigned r = 0; r < repeats; ++r) {
        std::fill(A.begin(), A.end(), 0.0); // важно: обнуляем A перед каждым прогоном
        Timer t; t.start();
        mul(A, B, C, n);
        acc += t.stop_sec();
    }
    return acc / repeats;
}

// -------------------------------
// Полиномиальная аппроксимация степени 3
// Возвращаем коэффициенты a3..a0 для P(n) = a3 n^3 + a2 n^2 + a1 n + a0
// Реализация через нормальные уравнения (быстро и достаточно для 5–10 точек)
// -------------------------------
struct Poly { double a3, a2, a1, a0; };

Poly polyfit3(const std::vector<double>& x, const std::vector<double>& y) {
    // Строим матрицу Вандермонда и решаем (X^T X) c = X^T y
    // Для краткости — минимальная реализация Гаусса 4×4 (без устойчивых библиотек).
    auto s = [&](int p){ double r=0; for(double xi: x) r+= std::pow(xi,p); return r; };
    double S0 = x.size();
    double S1 = s(1), S2 = s(2), S3 = s(3), S4 = s(4), S5 = s(5), S6 = s(6);
    double T0=0,T1=0,T2=0,T3=0;
    for (size_t i=0;i<x.size();++i) {
        T0 += y[i];
        T1 += y[i]*x[i];
        T2 += y[i]*x[i]*x[i];
        T3 += y[i]*x[i]*x[i]*x[i];
    }
    // Система:
    // [S6 S5 S4 S3][a3] = [T3]
    // [S5 S4 S3 S2][a2]   [T2]
    // [S4 S3 S2 S1][a1]   [T1]
    // [S3 S2 S1 S0][a0]   [T0]
    double A[4][5] = {
        {S6,S5,S4,S3,T3},
        {S5,S4,S3,S2,T2},
        {S4,S3,S2,S1,T1},
        {S3,S2,S1,S0,T0}
    };
    // Прямой ход Гаусса
    for(int i=0;i<4;i++){
        // pivot
        int piv=i;
        for(int r=i+1;r<4;r++) if (std::fabs(A[r][i])>std::fabs(A[piv][i])) piv=r;
        for(int c=i;c<5;c++) std::swap(A[i][c],A[piv][c]);
        // нормировка строки
        double d=A[i][i];
        for(int c=i;c<5;c++) A[i][c]/=d;
        // зануление ниже
        for(int r=i+1;r<4;r++){
            double k=A[r][i];
            for(int c=i;c<5;c++) A[r][c]-=k*A[i][c];
        }
    }
    // Обратный ход
    for(int i=3;i>=0;i--){
        for(int r=i-1;r>=0;r--){
            double k=A[r][i];
            A[r][i]=0;
            A[r][4]-=k*A[i][4];
        }
    }
    return Poly{A[0][4],A[1][4],A[2][4],A[3][4]};
}

// Коэффициент детерминации R^2 (качество аппроксимации)
double r2_score(const std::vector<double>& y, const std::vector<double>& yhat) {
    double mean=0; for(double v:y) mean+=v; mean/=y.size();
    double ss_tot=0, ss_res=0;
    for(size_t i=0;i<y.size();++i){ ss_tot+=(y[i]-mean)*(y[i]-mean); ss_res+=(y[i]-yhat[i])*(y[i]-yhat[i]); }
    return 1.0 - ss_res/ss_tot;
}

// Удобный форматтер под Excel RU (замена точки на запятую)
std::string fmt(double v) {
    std::ostringstream oss;
    oss.setf(std::ios::fixed); oss<<std::setprecision(6)<<v;
    auto s = oss.str();
    std::replace(s.begin(), s.end(), '.', ',');
    return s;
}

int main() {
    // Набор размеров матриц. Можно менять на свои.
    const std::vector<size_t> sizes = {80, 400, 700, 1000, 1500};
    // Повторы для усреднения (чем больше n — тем меньше повторов, чтобы не ждать долго)
    const unsigned repeats = 3;

    // Открываем CSV на запись (разделитель — ';')
    std::ofstream csv("results.csv");
    csv << "n;i,j,k;i,k,j;k,i,j;k,j,i;j,i,k;j,k,i\n";

    std::cout << "Matrix Multiplication Performance Comparison\n"
              << "==========================================\n\n";

    // Для сбора точек по каждому порядку (для последующей аппроксимации)
    std::vector<double> xs, y_ijk, y_ikj, y_kij, y_kji, y_jik, y_jki;

    for (size_t n : sizes) {
        Matrix A(n*n, 0.0), B(n*n), C(n*n);
        // Инициализация входных матриц одинаковым seed для сравнимости
        fill_random(B, n, /*seed*/12345);
        fill_random(C, n, /*seed*/67890);

        // Измеряем все порядки
        double t_ijk = measure_mul(mul_ijk, n, repeats, B, C, A);
        double t_ikj = measure_mul(mul_ikj, n, repeats, B, C, A);
        double t_kij = measure_mul(mul_kij, n, repeats, B, C, A);
        double t_kji = measure_mul(mul_kji, n, repeats, B, C, A);
        double t_jik = measure_mul(mul_jik, n, repeats, B, C, A);
        double t_jki = measure_mul(mul_jki, n, repeats, B, C, A);

        // Консольный отчёт по размеру
        std::cout << "Matrix size: " << n << "x" << n << "\n"
                  << "==========================================\n";
        std::cout << "(i,j,k): " << t_ijk << "s\n";
        std::cout << "(i,k,j): " << t_ikj << "s\n";
        std::cout << "(k,i,j): " << t_kij << "s\n";
        std::cout << "(k,j,i): " << t_kji << "s\n";
        std::cout << "(j,i,k): " << t_jik << "s\n";
        std::cout << "(j,k,i): " << t_jki << "s\n";

        // Поиск лучшего порядка для этого n
        double best = std::numeric_limits<double>::max();
        std::string best_name;
        auto upd_best = [&](double val, const char* name){
            if (val < best) { best = val; best_name = name; }
        };
        upd_best(t_ijk,"(i,k,j)"); // опечатка специально? нет — дальше правильно:
        best = std::numeric_limits<double>::max(); best_name.clear();
        upd_best(t_ijk,"(i,j,k)");
        upd_best(t_ikj,"(i,k,j)");
        upd_best(t_kij,"(k,i,j)");
        upd_best(t_kji,"(k,j,i)");
        upd_best(t_jik,"(j,i,k)");
        upd_best(t_jki,"(j,k,i)");
        std::cout << "Best method: " << best_name << " (" << best << "s)\n\n";

        // Пишем строку в CSV (десятичная запятая)
        csv << n << ';' << fmt(t_ijk) << ';' << fmt(t_ikj) << ';' << fmt(t_kij)
            << ';' << fmt(t_kji) << ';' << fmt(t_jik) << ';' << fmt(t_jki) << '\n';

        // Копим точки для трендов
        xs.push_back(static_cast<double>(n));
        y_ijk.push_back(t_ijk); y_ikj.push_back(t_ikj);
        y_kij.push_back(t_kij); y_kji.push_back(t_kji);
        y_jik.push_back(t_jik); y_jki.push_back(t_jki);
    }

    csv.close();
    std::cout << "Результаты сохранены в файл: results.csv\n\n";

    // Аппроксимация и печать формул + R^2 в консоль
    auto print_trend = [&](const char* name, const std::vector<double>& y) {
        Poly c = polyfit3(xs, y);
        std::vector<double> yhat; yhat.reserve(xs.size());
        for (double n : xs) yhat.push_back(c.a3*n*n*n + c.a2*n*n + c.a1*n + c.a0);
        double r2 = r2_score(y, yhat);

        std::cout.setf(std::ios::scientific);
        std::cout << name << " trend: y = "
                  << c.a3 << "*n^3 + " << c.a2 << "*n^2 + "
                  << c.a1 << "*n + " << c.a0
                  << " ; R^2 = " << r2 << "\n";
        std::cout.unsetf(std::ios::scientific);
    };

    std::cout << "=== Polynomial trend (3rd degree) for each order ===\n";
    print_trend("i,j,k", y_ijk);
    print_trend("i,k,j", y_ikj);
    print_trend("k,i,j", y_kij);
    print_trend("k,j,i", y_kji);
    print_trend("j,i,k", y_jik);
    print_trend("j,k,i", y_jki);

    // Средний тренд по всем порядкам (для демонстрации общего O(n^3))
    std::vector<double> y_avg(xs.size());
    for (size_t i=0;i<xs.size();++i)
        y_avg[i] = (y_ijk[i]+y_ikj[i]+y_kij[i]+y_kji[i]+y_jik[i]+y_jki[i]) / 6.0;

    std::cout << "\nОбщий тренд (среднее по всем порядкам):\n";
    print_trend("avg(all orders)", y_avg);

    std::cout << "\n=== Analysis (short) ===\n"
              << "- Перестановка циклов не меняет математику, но меняет паттерн доступа к памяти.\n"
              << "- Порядки, которые проходят по строкам/столбцам последовательно, лучше используют кэш.\n";
    return 0;
}
