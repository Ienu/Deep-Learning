#include <iostream>
#include "mat.h"

#define	MAX_ITR	7

using namespace std;

int main()
{
	Mat M_x(10, 3);
	M_x.data[0][0] = 1;
	M_x.data[1][0] = 1;
	M_x.data[2][0] = 1;
	M_x.data[3][0] = 1;
	M_x.data[4][0] = 1;
	M_x.data[5][0] = 1;
	M_x.data[6][0] = 1;
	M_x.data[7][0] = 1;
	M_x.data[8][0] = 1;
	M_x.data[9][0] = 1;

	M_x.data[0][1] = 20;	M_x.data[0][2] = 40;
	M_x.data[1][1] = 30;	M_x.data[1][2] = 50;
	M_x.data[2][1] = 20;	M_x.data[2][2] = 75;
	M_x.data[3][1] = 40;	M_x.data[3][2] = 65;
	M_x.data[4][1] = 40;	M_x.data[4][2] = 75;
	M_x.data[5][1] = 50;	M_x.data[5][2] = 80;
	M_x.data[6][1] = 50;	M_x.data[6][2] = 65;
	M_x.data[7][1] = 40;	M_x.data[7][2] = 70;
	M_x.data[8][1] = 30;	M_x.data[8][2] = 80;
	M_x.data[9][1] = 30;	M_x.data[9][2] = 65;

	Mat M_y(10, 1);
	M_y.data[0][0] = 0;
	M_y.data[1][0] = 0;
	M_y.data[2][0] = 0;
	M_y.data[3][0] = 0;
	M_y.data[4][0] = 0;
	M_y.data[5][0] = 1;
	M_y.data[6][0] = 1;
	M_y.data[7][0] = 1;
	M_y.data[8][0] = 1;
	M_y.data[9][0] = 1;

	Mat M_theta(3, 1);

	//Mat M_J(MAX_ITR, 1);

	for(int i = 0; i < MAX_ITR; i++)
	{
		Mat M_z(10, 1);
		M_z = M_x * M_theta;

		Mat M_h(10, 1);
		M_h = MSigmoid(M_z);

		Mat M_grad(3, 1);
		M_grad = (1.0 / 10.0) * MTrans(M_x) * (M_h - M_y);

		Mat M_H(3, 3);
		M_H = (1.0 / 10.0) * MTrans(M_x) * MDiag(M_h) * MDiag(1 - M_h) * M_x;
		//M_J.data[i][0] = (1.0 / 10.0) *  
		M_theta = M_theta - MInv(M_H) * M_grad;
	}

	M_theta.MShow();

	return 0;
}
