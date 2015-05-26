#include <iostream>
#include <cv.h>
#include <highgui.h>

#include "mat.h"

#pragma comment(lib, "cv.lib")
#pragma comment(lib, "highgui.lib")
#pragma comment(lib, "cxcore.lib")



using namespace std;

int main()
{
	double** pca_data = MReadFile("pcaData.txt", 2, 45);
	Mat M_data(2, 45, pca_data);
	M_data.MShow();

	Mat M_sigma(2, 2);
	M_sigma = (1.0 / 45.0) * M_data * MTrans(M_data);
	M_sigma.MShow();

	double d[4];
	d[0] = M_sigma[0][0];
	d[1] = M_sigma[0][1];
	d[2] = M_sigma[1][0];
	d[3] = M_sigma[1][1];

	CvMat* A = cvCreateMat(2, 2, CV_64F);
	CvMat* W = cvCreateMat(2, 2, CV_64F);
	CvMat* U = cvCreateMat(2, 2, CV_64F);
	CvMat* V = cvCreateMat(2, 2, CV_64F);

	Mat M_U(2, 2);
	Mat M_S(2, 2);

	for(int i = 0; i < 4; i++)
	{
		A->data.db[i] = d[i];
	}

	cvSVD(A, W, U, V);

	cout << "U" << endl;
	for(i = 0; i < 4; i++)
	{
		cout << U->data.db[i] << endl;
	}
	M_U[0][0] = U->data.db[0];
	M_U[0][1] = U->data.db[1];
	M_U[1][0] = U->data.db[2];
	M_U[1][1] = U->data.db[3];

	M_S[0][0] = W->data.db[0];
	M_S[0][1] = W->data.db[1];
	M_S[1][0] = W->data.db[2];
	M_S[1][1] = W->data.db[3];

	cout << "V" << endl;
	for(i = 0; i < 4; i++)
	{
		cout << V->data.db[i] << endl;
	}

	cout << "S" << endl;
	for(i = 0; i < 4; i++)
	{
		cout << W->data.db[i] << endl;
	}

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
