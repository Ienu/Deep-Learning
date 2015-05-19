#include <iostream>
#include <fstream>
#include <cmath>
#include "mat.h"

#define MAX_ITR	15

using namespace std;

int main()
{
	double** d_x = MReadFile("ex5Logx.txt", 117, 2);
	double** x = new double*[117];
	for(int i = 0; i < 117; i++)
	{
		x[i] = new double[28];
		for(int j = 0; j < 28; j++)
		{
			x[i][j] =  1;
		}
	}
	for(int k = 0; k < 117; k++)
	{
		int index = 1;
		for(i = 1; i <= 6; i++)
		{
			for(int j = 0; j <= i; j++)
			{
				x[k][index] = pow(d_x[k][0], i - j) * pow(d_x[k][1], j);
				index++;
			}
		}
	}
	Mat M_x(117, 28, x);
	double** d_y = MReadFile("ex5Logy.txt", 117, 1);
	Mat M_y(117, 1, d_y);

	Mat M_theta(28, 1);

	double lambda = 1.0;

	for(i = 0; i < MAX_ITR; i++)
	{
		Mat M_z(117, 1);
		M_z = M_x * M_theta;

		Mat M_h(117, 1);
		M_h = MSigmoid(M_z);		

		Mat M_G(28, 1);
		M_G = (lambda / 117.0) * M_theta;
		M_G[0][0] = 0;

		Mat M_L(28, 28);
		M_L = (lambda / 117.0) * MEye(28);
		M_L[0][0] = 0;

		Mat M_grad(28, 1);
		M_grad = ((1.0 / 117.0) * MTrans(M_x) * (M_h - M_y)) + M_G;

		Mat M_H(28, 28);
		M_H = ((1.0 / 117.0) * MTrans(M_x) * MDiag(M_h) * MDiag(1 - M_h) * M_x) + M_L;

		M_theta = M_theta - MInv(M_H) * M_grad;
	}

	M_theta.MShow();

	return 0;
}
