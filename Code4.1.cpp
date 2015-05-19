#include <iostream>
#include <cmath>
#include "mat.h"

using namespace std;

int main()
{
	Mat M_x(7, 1);
	M_x.data[0][0] = -1;
	M_x.data[1][0] = -0.7;
	M_x.data[2][0] = -0.4;
	M_x.data[3][0] = -0.1;
	M_x.data[4][0] = 0.2;
	M_x.data[5][0] = 0.5;
	M_x.data[6][0] = 0.8;

	Mat M_y(7, 1);
	M_y.data[0][0] = 2.1;
	M_y.data[1][0] = 1.2;
	M_y.data[2][0] = 0.4;
	M_y.data[3][0] = 0.5;
	M_y.data[4][0] = 0.5;
	M_y.data[5][0] = 0.2;
	M_y.data[6][0] = -0.3;

	Mat M_sx(7, 6);
	for(int i = 0; i < 7; i++)
	{
		for(int j = 0; j < 6; j++)
		{
			M_sx.data[i][j] = pow(M_x.data[i][0], j);
		}
	}

	Mat M_rm(6, 6);
	for(i = 1; i < 6; i++)
	{
		M_rm.data[i][i] = 1;
	}

	double d_lambda[3] = {0, 1, 10};

	Mat** pM_sida = new Mat*[3];
	for(i = 0; i < 3; i++)
	{
		pM_sida[i] = new Mat(6, 1);
	}

	for(i = 0; i < 3; i++)
	{
		*(pM_sida[i]) = MInv(MTrans(M_sx) * M_sx + d_lambda[i] * M_rm) * MTrans(M_sx) * M_y;
		pM_sida[i]->MShow();
	}

	delete[] pM_sida;

	return 0;
}
