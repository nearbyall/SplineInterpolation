#include "SplineInterpolation.h"

double** SplineInterpolation::readMatrixFromFile(std::string matrixpath, int n, int m)
{
    double** matrix = new double* [n];
    for (int i = 0; i < n; i++)
        matrix[i] = new double[m];
    std::ifstream in("interpolation_polynomial_41.lb3");

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            in >> matrix[i][j];
        }
    }

    return matrix;
}

void SplineInterpolation::printValues(double* x, double* y, int n)
{
    std::cout << "X:" << std::endl;
    for (int i = 0; i < n; i++)
        std::cout << x[i] << std::endl;
    std::cout << "Y:" << std::endl;
    for (int i = 0; i < n; i++)
        std::cout << y[i] << std::endl;
}

double SplineInterpolation::function(double x2, double k1)
{
    double y;
    y = 1. / (1. + k1 * x2 * x2);
    return y;
}

void SplineInterpolation::sweepMethod(double* x, double* y, double* h, double* w, double* u, double* v, double* k, double* l, double* c, int n, int m)
{
    for (int i = 1; i <= n; i++) {
        h[i] = x[i] - x[i - 1];
    }

    for (int i = 2; i <= n; i++) {
        w[i] = h[i - 1];
        u[i] = 2 * (h[i - 1] + h[i]);
        v[i] = 3 * (((y[i] - y[i - 1]) / h[i]) - ((y[i - 1] - y[i - 2]) / h[i - 1]));
    }

    k[1] = 0;
    l[1] = 0;

    for (int i = 2; i <= n; i++) {
        k[i] = (v[i] - w[i] * k[i - 1]) / (v[i] - w[i] * k[i - 1]);
        l[i] = h[i] / (v[i] - w[i] * k[i - 1]);
    }

    c[n + 1] = 0;
    for (int i = n; i >= 1; i--) {
        c[i] = k[i] - l[i] * c[i + 1];
    }

    for (int i = 1; i <= n; i++) {
        std::cout << c[i] << std::endl;
    }
}
