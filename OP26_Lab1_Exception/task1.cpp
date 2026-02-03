#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <stdexcept>
#include <algorithm>

using namespace std;

// === Винятки для перемикання між алгоритмами ===
class SwitchToAlg2 : public exception {};
class SwitchToAlg3 : public exception {};
class SwitchToAlg4 : public exception {}; // Перерахунок через Алг 4
class SwitchToAlg5 : public exception {}; // Глобальна помилка

// === Структура даних для таблиць ===
struct Point {
    double x, t, u;
};

// Дані Таблиці 1 (dat_X_1_1.dat) - для |x| <= 1 [cite: 336]
const vector<Point> table1 = {
    {-1.000,-4.935,1.935}, {-0.900,-3.013,0.464}, {-0.800,-2.316,1.327},
    {-0.700,-1.819,1.976}, {-0.600,-1.425,2.502}, {-0.500,-1.097,2.951},
    {-0.400,-0.816,3.344}, {-0.300,-0.571,3.695}, {-0.200,-0.357,4.013},
    {-0.100,-0.167,4.303}, {0.000,0.000,4.571},   {0.100,0.147,4.618},
    {0.200,0.276,4.645},   {0.300,0.386,4.652},   {0.400,0.477,4.636},
    {0.500,0.548,4.596},   {0.600,0.597,4.524},   {0.700,0.617,4.412},
    {0.800,0.597,4.240},   {0.900,0.505,3.956},   {1.000,0.000,3.000}
};

// Дані Таблиці 2 (dat_X_1_00.dat) - для x < -1 (x = -1/x) [cite: 337] (адаптовано за логікою джерел)
const vector<Point> table2 = {
    {0.000,-4.935,1.935}, {0.050,-2.663,1.885}, {0.100,-1.618,1.834},
    {0.150,-0.773,1.784}, {0.200,-0.034,1.732}, {0.250,0.635,1.679},
    {0.300,1.253,1.625},  {0.350,1.829,1.570},  {0.400,2.369,1.512},
    {0.450,2.877,1.452},  {0.500,3.356,1.388},  {0.550,3.806,1.322},
    {0.600,4.228,1.251},  {0.650,4.622,1.175},  {0.700,4.987,1.093},
    {0.750,5.320,1.003},  {0.800,5.618,0.905},  {0.850,5.876,0.796},
    {0.900,6.080,0.675},  {0.950,6.199,0.536},  {1.000,5.890,0.377}
};

// Дані Таблиці 3 (dat_X_00_1.dat) - для x > 1 (x = -1/x) [cite: 335]
const vector<Point> table3 = {
    {0.000,-4.935,1.935},  {-0.050,-4.435,1.835}, {-0.100,-3.936,1.735},
    {-0.150,-3.440,1.636}, {-0.200,1.537,-2.948}, {-0.250,-2.461,1.440},
    {-0.300,-1.980,1.344}, {-0.350,-1.506,1.249}, {-0.400,-1.041,1.156},
    {-0.450,-0.585,1.065}, {-0.500,0.976,-0.141}, {-0.550,0.292,0.889},
    {-0.600,0.712,0.806},  {-0.650,1.117,0.724},  {-0.700,1.507,0.646},
    {-0.750,0.572,1.882},  {-0.800,2.239,0.500},  {-0.850,2.578,0.432},
    {-0.900,2.898,0.368},  {-0.950,3.199,0.308},  {-1.000,3.480,0.252}
};

// === Допоміжні функції ===

// Інтерполяція [cite: 415, 416]
pair<double, double> interpolate(double x, const vector<Point>& table) {
    if (table.empty()) throw SwitchToAlg5();
    
    // Проста перевірка меж та пошук інтервалу
    size_t n = table.size();
    
    // Якщо масив відсортований за спаданням (Таблиця 3) або зростанням
    bool ascending = table[0].x < table[n-1].x;

    for (size_t i = 0; i < n - 1; ++i) {
        double x1 = table[i].x;
        double x2 = table[i+1].x;
        
        bool in_range = ascending ? (x >= x1 && x <= x2) : (x <= x1 && x >= x2);

        if (in_range) {
            double t1 = table[i].t, t2 = table[i+1].t;
            double u1 = table[i].u, u2 = table[i+1].u;
            
            if (abs(x2 - x1) < 1e-9) return {t1, u1};
            
            double t_val = t1 + (t2 - t1) * (x - x1) / (x2 - x1);
            double u_val = u1 + (u2 - u1) * (x - x1) / (x2 - x1);
            return {t_val, u_val};
        }
    }
    // Повертаємо найближче граничне значення, якщо вийшли за межі (або можна кидати Alg5)
    if (ascending) return (x < table[0].x) ? make_pair(table[0].t, table[0].u) : make_pair(table[n-1].t, table[n-1].u);
    else return (x > table[0].x) ? make_pair(table[0].t, table[0].u) : make_pair(table[n-1].t, table[n-1].u);
}

// Отримання T(x) та U(x) [cite: 408-411]
pair<double, double> get_TU(double x) {
    if (abs(x) <= 1.0) return interpolate(x, table1);
    else if (x < -1.0) return interpolate(-1.0 / x, table2);
    else return interpolate(-1.0 / x, table3);
}

double T(double x) { return get_TU(x).first; }
double U(double x) { return get_TU(x).second; }

double Srz(double x, double y, double z) {
    // [cite: 407]
    if (x > y) return T(x) + U(z) - T(y);
    else return T(y) + U(y) - U(z);
}

// === Прототипи функцій ===
double Rrz2_Alg4(double x, double y, double z);
double Krm_Alg4(double x, double y, double z);
double Rrz_Alg2(double x, double y, double z);
double Rrz_Alg3(double x, double y, double z);

// === АЛГОРИТМ 1 (Основний) ===

// [cite: 402, 403]
double Srs_Alg1(double x, double y, double z) {
    if (z > y) {
        double val = z*z + x*y;
        if (val > 0) return Srz(x, y, z) + y * sqrt(val);
        else throw SwitchToAlg2(); // Пункт 5.1
    } else {
        double val = x*x + z*y;
        if (val > 0) return y + Srz(z, x, y) * sqrt(val);
        else throw SwitchToAlg3(); // Пункт 5.2
    }
}

// [cite: 404, 406]
double Srs1_Alg1(double x, double y, double z) {
    if (z > y) {
        double val = z*z + x*y;
        if (val > 1) return Srz(x, y, z) + y * log(val);
        else throw SwitchToAlg2(); // Пункт 6.1
    } else {
        double val = x*x + z*y;
        if (val > 1) return y + Srz(z, x, y) * sqrt(val);
        else throw SwitchToAlg4(); // Пункт 6.2 (Перерахунок через Алг 4)
    }
}

// [cite: 401]
double Qrz_Alg1(double x, double y) {
    if (abs(x) < 1) return x * Srs_Alg1(x, y, x);
    else return y * Srs1_Alg1(y, x, y);
}

// [cite: 400]
double Rrz_Alg1(double x, double y, double z) {
    try {
        if (x > y) return x * z * Qrz_Alg1(y, z) - x;
        else return y * x * Qrz_Alg1(x, y) + y;
    } catch (SwitchToAlg2&) { return Rrz_Alg2(x, y, z); }
      catch (SwitchToAlg3&) { return Rrz_Alg3(x, y, z); }
      catch (SwitchToAlg4&) { return Rrz2_Alg4(x, y, z); } // Використовуємо Rrz2 з Алг 4
}

// 
double Km_Alg1(double x, double y, double z) {
    return 73.1389 * Rrz_Alg1(x, y, y) + 14.838 * Rrz_Alg1(x - y, z, y);
}

// === АЛГОРИТМ 2 [cite: 417-420] ===
double Srs1_Alg2(double x, double y, double z) {
    if (z > y) return Srz(x, y, z) + 1.44 * y * z;
    else return y + 1.44 * Srz(z, x, y);
}
double Qrz1_Alg2(double x, double y) {
    if (abs(y) < 1) return x * Srs1_Alg2(x, y, x);
    else return y * Srs1_Alg2(y, x, y);
}
double Rrz_Alg2(double x, double y, double z) {
    if (x > y) return x * y * Qrz1_Alg2(y, z);
    else return x * z * Qrz1_Alg2(x, y);
}

// === АЛГОРИТМ 3 [cite: 422-424] ===
double Srs2_Alg3(double x, double y, double z) {
    if (z > y) return Srz(x, y, z) + y * x;
    else return y * z + Srz(z, x, y);
}
double Qrz2_Alg3(double x, double y) {
    if (abs(x) < 1) return x * Srs2_Alg3(x, y, x);
    else return y * Srs2_Alg3(y, x, y);
}
double Rrz_Alg3(double x, double y, double z) {
    if (x > y) return x * y * Qrz2_Alg3(y, z);
    else return y * z * Qrz2_Alg3(x, y);
}

// === АЛГОРИТМ 4 [cite: 425-429] ===
double Srs3_Alg4(double x, double y, double z) {
    if (z > y) return Srz(x, y, z) + y * x;
    else return y * z + Srz(z, x, y);
}
double Qrz3_Alg4(double x, double y) {
    if (abs(x) < 1) return x * Srs3_Alg4(x, y, x);
    else return x * Srs3_Alg4(y, x, y);
}
double Rrz2_Alg4(double x, double y, double z) {
    if (x > y) return y * Qrz3_Alg4(y, z);
    else return z * Qrz3_Alg4(x, y);
}
double Krm_Alg4(double x, double y, double z) {
    return 83.1389 * Rrz2_Alg4(x, y, z) + 4.838 * Rrz2_Alg4(x, z, y);
}

// === АЛГОРИТМ 5 (Fallback) [cite: 431] ===
double fun_Alg5(double x, double y, double z) {
    return 4.349 * x * z + 23.23 * y - 2.348 * x * y * z;
}

// === Головна функція розрахунку ===
double calculate_fun(double x, double y, double z) {
    try {
        //  fun = x*Km + y*Krm - z*Krm
        // Km з Алг 1, Krm беремо з Алг 4 (бо в Алг 1 він не визначений, але використовується)
        double km = Km_Alg1(x, y, z);
        double krm = Krm_Alg4(x, z, y); // Увага на порядок аргументів у виклику з Алг 1
        
        return x * km + y * krm - z * krm;
    } catch (SwitchToAlg5&) {
        return fun_Alg5(x, y, z); // Якщо проблеми з таблицями
    } catch (...) {
        return fun_Alg5(x, y, z); // Інші помилки
    }
}

int main() {
    double x, y, z;
    cout << "Vvedit x, y, z: ";
    if (cin >> x >> y >> z) {
        double result = calculate_fun(x, y, z);
        cout << "Result fun(x, y, z) = " << result << endl;
    } else {
        cout << "Nekorektni dani." << endl;
    }
    return 0;
}