#include <iostream>
#include "Mat.h"

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

	Mat M_w(2, 1);
	M_w = MInv(MTrans(M_x) * M_x) * MTrans(M_x) * M_y;

	M_w.MShow();

	return 0;
}
