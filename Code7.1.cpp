#include <iostream>
#include "mat.h"

using namespace std;

int main()
{
	double** pca_data = MReadFile("pcaData.txt", 2, 45);
	Mat M_data(2, 45, pca_data);
	M_data.MShow();

	Mat M_sigma(2, 2);
	M_sigma = (1.0 / 45.0) * M_data * MTrans(M_data);
	M_sigma.MShow();

	Mat M_U(2, 2);
	Mat M_V(2, 2);
	Mat M_W(2, 2);
	Mat M_S(2, 2);

	M_S = M_sigma.SVD(M_U, M_W, M_V);

	Mat M_xRot(2, 45);
	M_xRot = MTrans(M_U) * M_data;
	M_xRot.MShow();

	double epsilon = 1e-5;

	Mat M_xPCAWhite(2, 45);
	M_xPCAWhite = MDiag(1.0 / MDPow(M_S * MOnes(2, 1) + epsilon, 0.5)) * MTrans(M_U) * M_data;
	M_xPCAWhite.MShow();

	Mat M_xZCAWhite(2, 45);
	M_xZCAWhite = M_U * M_xPCAWhite;
	M_xZCAWhite.MShow();

	return 0;
}
