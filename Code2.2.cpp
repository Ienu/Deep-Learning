#include <iostream>
#include "mat.h"

#define MAX_ITR	1500

using namespace std;

int main()
{
	Mat M_x(5, 2);
	M_x.data[0][0] = 1;
	M_x.data[1][0] = 1;
	M_x.data[2][0] = 1;
	M_x.data[3][0] = 1;
	M_x.data[4][0] = 1;
	
	M_x.data[0][1] = 2;
	M_x.data[1][1] = 3;
	M_x.data[2][1] = 4;
	M_x.data[3][1] = 5;
	M_x.data[4][1] = 6;
	
	Mat M_y(5, 1);
	M_y.data[0][0] = 0.82;
	M_y.data[1][0] = 0.93;
	M_y.data[2][0] = 1.24;
	M_y.data[3][0] = 1.22;
	M_y.data[4][0] = 1.45;

	Mat M_theta(2, 1);

	double d_alpha = 0.07;

	for(int i = 0; i < MAX_ITR; i++)
	{
		Mat M_grad(2, 1);
		M_grad = (1.0 / 5.0) * MTrans(M_x) * (M_x * M_theta - M_y);
		M_theta = M_theta - d_alpha * M_grad;
	}

	M_theta.MShow();

	return 0;
}
