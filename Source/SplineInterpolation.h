#pragma once

#include <iostream>
#include <string>
#include <fstream>

class SplineInterpolation
{

private:

	SplineInterpolation() {}

public:

	//Считывание матрицы из файла
	static double** readMatrixFromFile(std::string matrixpath, int n, int m);

	//Вывод значений x и y в консоль
	static void printValues(double* x, double* y, int n);

	//Выбранная функция
	static double function(double x2, double k1);

	//Метод прогонки
	static void sweepMethod(double* x, double* y, double* h, double* w, double* u, double* v, double* k, double* l, double* c, int n, int m);

};

