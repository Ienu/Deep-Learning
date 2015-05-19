// MAT_H
#ifndef MAT_H
#define MAT_H

#include <iostream>
#include <cmath>

using namespace std;

class Mat
{
public:
	double** data;
	
public:
	int row;
	int col;
	
public:
	
	Mat& operator = (const Mat& M);	
	double* & operator [] (const int& index);
	void MShow();
	bool QR(Mat& Q, Mat& R);
	
	Mat(int row, int col);
	~Mat();
	
};
Mat& MEye(const int& n)
{
	if(n <= 0)
	{
		cerr << "n <= 0" << endl;
		exit(0);
	}
	Mat* pM_res = new Mat(n, n);
	for(int i = 0; i < n; i++)
	{
		pM_res->data[i][i] = 1;
	}
	return *pM_res;
}

double* & Mat::operator [] (const int& index)
{
	if(data == NULL)
	{
		cerr << "Matrix has no data" << endl;
		exit(0);
	}
	if(index < 0 || index >= row)
	{
		cerr << "Index exceed dimension" << endl;
		exit(0);
	}
	return data[index];
}
Mat::Mat(int row, int col)
{
	if(row <= 0 || col <= 0)
	{
		cerr << "row or col <= 0" << endl;
		exit(0);
	}
	this->row = row;
	this->col = col;
	
	data = new double*[row];
	for(int i = 0; i < row; i++)
	{
		data[i] = new double[col];
		for(int j = 0; j < col; j++)
		{
			data[i][j] = 0;
		}
	}
}
Mat::~Mat()
{
	if(data == NULL)
	{
		exit(0);
	}
	for(int i = 0; i < row; i++)
	{
		if(data[i] == NULL)
		{
			exit(0);
		}
		delete[] data[i];
	}
	delete[] data;
}
void Mat::MShow()
{
	for(int i = 0; i < row; i++)
	{
		for(int j = 0; j < col; j++)
		{
			cout << data[i][j] << "\t";
		}
		cout << endl;
	}
}
Mat& Mat::operator =(const Mat& M)
{
	if(M.data == NULL)
	{
		cerr << "Matrix has no data" << endl;
		exit(0);
	}
	if(M.row != row || M.col != col)
	{
		cerr << "Matrix dimension do not match" << endl;
		exit(0);
	}
	for(int i = 0; i < row; i++)
	{
		for(int j = 0; j < col; j++)
		{
			data[i][j] = M.data[i][j];
		}
	}
	return *this;
}
Mat& MInv(const Mat& M)
{
	if(M.data == NULL)
	{
		cerr << "Matrix has no data" << endl;
		exit(0);
	}
	if(M.col != M.row)
	{
		cerr << "Matrix is not square" << endl;
		exit(0);
	}
	double d_value;
	int i_n = M.row;
	
	Mat* pM_res = new Mat(i_n, i_n);
	*pM_res = M;
	
	for(int i = 0; i < i_n; i++)
	{
		d_value = 1.0 / pM_res->data[i][i];
		pM_res->data[i][i] = d_value;
		for(int j = 0; j < i_n; j++)
		{
			if(i != j)
			{
				pM_res->data[i][j] *= -d_value;
			}
		}
		for(j = 0; j < i_n; j++)
		{
			if(i != j)
			{
				pM_res->data[j][i] *= d_value;
			}
		}
		for(j = 0; j < i_n; j++)
		{
			if(i != j)
			{
				for(int k = 0; k < i_n; k++)
				{
					if(i != k)
					{
						pM_res->data[j][k] += pM_res->data[j][i] * pM_res->data[i][k] / d_value;
					}
				}
			}
		}
	}
	return *pM_res;
}

Mat& operator +(const Mat& M1, const Mat& M2)
{
	if(M1.data == NULL || M2.data == NULL)
	{
		cerr << "M1 or M2 has no data" << endl;
		exit(0);
	}
	if(M1.row != M2.row || M1.col != M2.col)
	{
		cerr << "Matrix dimension do not match" << endl;
		exit(0);
	}
	Mat* pM_res = new Mat(M1.row, M1.col);
	for(int i = 0; i < pM_res->row; i++)
	{
		for(int j = 0; j < pM_res->col; j++)
		{
			pM_res->data[i][j] = M1.data[i][j] + M2.data[i][j];
		}
	}
	return *pM_res;
}

Mat& operator -(const Mat& M1, const Mat& M2)
{
	if(M1.data == NULL || M2.data == NULL)
	{
		cerr << "M1 or M2 has no data" << endl;
		exit(0);
	}
	if(M1.row != M2.row || M1.col != M2.col)
	{
		cerr << "Matrix dimension do not match" << endl;
		exit(0);
	}
	Mat* pM_res = new Mat(M1.row, M1.col);
	for(int i = 0; i < pM_res->row; i++)
	{
		for(int j = 0; j < pM_res->col; j++)
		{
			pM_res->data[i][j] = M1.data[i][j] - M2.data[i][j];
		}
	}
	return *pM_res;
}

Mat& operator -(const double& s, const Mat& M2)
{
	if(M2.data == NULL)
	{
		cerr << "M2 has no data" << endl;
		exit(0);
	}
	Mat* pM_res = new Mat(M2.row, M2.col);
	for(int i = 0; i < pM_res->row; i++)
	{
		for(int j = 0; j < pM_res->col; j++)
		{
			pM_res->data[i][j] = s - M2.data[i][j];
		}
	}
	return *pM_res;
}

Mat& operator *(const Mat& M1, const Mat& M2)
{
	if(M1.data == NULL || M2.data == NULL)
	{
		cerr << "M1 or M2 has no data" << endl;
		exit(0);
	}
	if(M1.col != M2.row)
	{
		cerr << "Matrix dimension do not match" << endl;
		exit(0);
	}
	Mat* pM_res = new Mat(M1.row, M2.col);
	for(int i = 0; i < pM_res->row; i++)
	{
		for(int j = 0; j < pM_res->col; j++)
		{
			for(int k = 0; k < M1.col; k++)
			{
				pM_res->data[i][j] += M1.data[i][k] * M2.data[k][j];
			}
		}
	}
	return *pM_res;
}

Mat& operator *(const double& s, const Mat& M2)
{
	if(M2.data == NULL)
	{
		cerr << "M2 has no data" << endl;
		exit(0);
	}
	
	Mat* pM_res = new Mat(M2.row, M2.col);
	for(int i = 0; i < pM_res->row; i++)
	{
		for(int j = 0; j < pM_res->col; j++)
		{
			pM_res->data[i][j] = s * M2.data[i][j];
		}
	}
	return *pM_res;
}
Mat& MTrans(const Mat& M)
{
	if(M.data == NULL)
	{
		cerr << "Matrix has no data" << endl;
		exit(0);
	}
	Mat *pM_res = new Mat(M.col, M.row);
	for(int i = 0; i < pM_res->row; i++)
	{
		for(int j = 0; j < pM_res->col; j++)
		{
			pM_res->data[i][j] = M.data[j][i];
		}
	}
	return *pM_res;
}

Mat& MSigmoid(const Mat& M)
{
	if(M.data == NULL)
	{
		cerr << "Matrix has no data" << endl;
		exit(0);
	}
	Mat *pM_res = new Mat(M.row, M.col);
	for(int i = 0; i < pM_res->row; i++)
	{
		for(int j = 0; j < pM_res->col; j++)
		{
			pM_res->data[i][j] = 1.0 / (1.0 + exp(-M.data[i][j]));
		}
	}
	return *pM_res;
}

Mat& MDiag(const Mat& M)
{
	if(M.data == NULL)
	{
		cerr << "Matrix has no data" << endl;
		exit(0);
	}
	if(M.col != 1)
	{
		cerr << "Matrix is not cvector" << endl;
		exit(0);
	}
	Mat* pMat_res = new Mat(M.row, M.row);
	for(int i = 0; i < pMat_res->row; i++)
	{
		pMat_res->data[i][i] = M.data[i][0];
	}
	return *pMat_res;
}

bool Mat::QR(Mat& Q, Mat& R)
{
	if(data == NULL)
	{
		cerr << "Matrix has no data" << endl;
		exit(0);
	}
	if(Q.data == NULL || R.data == NULL)
	{
		cerr << "Q or R has not init data" << endl;
		exit(0);
	}
	if(row != col)
	{
		cerr << "Matrix is not square" << endl;
		exit(0);
	}
	if(Q.row != Q.col || R.row != R.col || Q.row != row)
	{
		cerr << "Matrix dimension do not match" << endl;
		exit(0);
	}
	
	int n = row;
	Mat M_A(n, n);
	M_A = *this;
	Q = MEye(n);
	Mat M_H(n, n);

	for(int k = 0; k < n; k++)
	{
		double sigma = -1;
		if(M_A[k][k] >= 0)
		{
			sigma = 1;
		}
		double sum = 0;
		for(int i = k; i < n; i++)
		{
			sum += M_A[i][k] * M_A[i][k];
		}
		sigma *= sqrt(sum);
		double rou = sigma * (sigma + M_A[k][k]);
		Mat M_U(n, 1);
		M_U[k][0] = sigma + M_A[k][k];
		for(int j = k + 1; j < n; j++)
		{
			M_U[j][0] = M_A[j][k];
		}
		M_H = MEye(n) - 1.0 / rou * M_U * MTrans(M_U);
		M_A = M_H * M_A;
		Q = Q * M_H;
	}
	R = M_A;
	
	return true;
}

#endif
