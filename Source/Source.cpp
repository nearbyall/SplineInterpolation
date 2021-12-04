#include <math.h>

#include "SplineInterpolation.h"

int main()
{

    const int n = 11;
    const int m = 2;
    double* S = new double[n];

    double* w = new double[n + 1];
    double* u = new double[n + 1];
    double* v = new double[n + 1];
    double* h = new double[n + 1];
    double* c = new double[n + 2];
    double* k = new double[n + 1];
    double* l = new double[n + 1];
    double* a = new double[n + 1];
    double* b = new double[n + 1];
    double* d = new double[n + 1];

    double* x = new double[n];
    double* y = new double[n];

    double** matrix = SplineInterpolation::readMatrixFromFile("interpolation_polynomial_41.lb3", n, m);

    for (int i = 0; i < n; i++) {
        x[i] = matrix[i][0];
        y[i] = matrix[i][1];
    }

    SplineInterpolation::printValues(x, y, n);

    // Метод прогонки и накождение коэффицентов a, b, c, d

    SplineInterpolation::sweepMethod(x, y, h, w, u, v, k, l, c, n, m);

    std::cout << std::endl;

    for (double z = -1.9; z < 1.9; z = z + 0.0376) {
        for (int i = 1; i <= n; i++) {
            if (z >= x[i - 1] && z < x[i]) {
                a[i] = y[i - 1];
                d[i] = (c[i + 1] + c[i]) / (3 * h[i]);
                b[i] = (y[i] - y[i - 1]) / h[i] - (h[i] * (c[i + 1] + 2 * c[i])) / 3;
                S[i] = a[i] + b[i] * (z - x[i - 1]) + c[i] * (z - x[i - 1]) * (z - x[i - 1]) + b[i] * (z - x[i - 1]) * (z - x[i - 1]) * (z - x[i - 1]);
                std::cout << z << "\t" << S[i] << std::endl;
            }
        }
    }

    std::ofstream out("Lab4_spline_Melnikov_41.cmp ");
    for (int i = 1; i < n; i++) {
        out << "a[" << i << "]: " << a[i] << "\t";
        out << "b[" << i << "]: " << b[i] << "\t";
        out << "c[" << i << "]: " << c[i] << "\t";
        out << "d[" << i << "]: " << d[i] << std::endl;
    }
    out.close();

    // Исследование сходимости интерполяционного процесса для функции 1/(1+kx^2)
    double* x2 = new double[n];
    double* y2 = new double[n];

    double k1 = 2;

    for (int i = 0; i < n; i++) {
        x2[i] = x[i];
        y2[i] = SplineInterpolation::function(x2[i], k1);
    }

    std::cout << std::endl;
    SplineInterpolation::printValues(x2, y2, n);
    std::cout << std::endl;

    SplineInterpolation::sweepMethod(x2, y2, h, w, u, v, k, l, c, n, m);

    // Вывод иксов, игриков, посчитанных с помощью интерполяции и игриков, посчитанных самой функцией.
    // Расчет максимальной ошибки
    out.open("Lab4_max_error_Melnikov_41.cmp");
    double max_error = abs(S[1] - SplineInterpolation::function(-2, k1));
    for (double z = -2; z < 2.04; z = z + 0.04) {
        for (int i = 1; i <= n; i++) {
            if (z >= x2[i - 1] && z < x2[i]) {
                a[i] = y2[i - 1];
                d[i] = (c[i + 1] + c[i]) / (3 * h[i]);
                b[i] = (y2[i] - y2[i - 1]) / h[i] - (h[i] * (c[i + 1] + 2 * c[i])) / 3;
                S[i] = a[i] + b[i] * (z - x2[i - 1]) + c[i] * (z - x2[i - 1]) * (z - x2[i - 1]) + b[i] * (z - x2[i - 1]) * (z - x2[i - 1]) * (z - x2[i - 1]);
                out << "X: " << z << "\t" << "\t" << "Y: (S) " << S[i] << "\t" << "\t";
                if (abs(S[i] - SplineInterpolation::function(z, k1)) > max_error) {
                    max_error = abs(S[i] - SplineInterpolation::function(z, k1));
                }
            }
        }
        out << "Y: " << SplineInterpolation::function(z, k1) << std::endl;
    }
    out << std::endl << "Max error = " << max_error << std::endl;
    out.close();

    return 0;
}