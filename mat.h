// MAT_H
#ifndef MAT_H
#define MAT_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "svd94c.h"

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
	Mat& operator *= (const double& s);
	Mat& operator /= (const Mat& M);
	Mat& operator += (const Mat& M);
	Mat& operator -= (const Mat& M);

	Mat All();
	
	double* & operator [] (const int& index);
	
	Mat operator - ();
	void MShow();
	bool QR(Mat& Q, Mat& R);
	Mat SVD(Mat& U, Mat& W, Mat& V);
	Mat Diag();
	
	Mat(const int& row, const int& col);
	Mat(const int& row, const int& col, double** d_data);

	Mat(const Mat& M);
	//Mat(double& d);
	
	Mat MSub(const int& li, const int& lj, const int& srow, const int& scol);
	void MPub(const Mat& M, const int& li, const int& lj);
	
	~Mat();
	
};

Mat Mat::All()
{
	if(data == NULL)
	{
		cerr << "All(): Matrix has no data" << endl;
		exit(0);
	}
	Mat M_res(row * col, 1);
	for(int i = 0; i < col; i++)
	{
		for(int j = 0; j < row; j++)
		{
			M_res.data[i * row + j][0] = data[j][i];
		}
	}
	return M_res;
}

Mat operator == (const Mat& M, const double& s)
{
	if(M.data == NULL)
	{
		cerr << "operator == (): Matrix has no data" << endl;
		exit(0);
	}

	Mat M_res(M.row, M.col);
	for(int i = 0; i < M.row; i++)
	{
		for(int j = 0; j < M.col; j++)
		{
			M_res.data[i][j] = (M.data[i][j] == s);
		}
	}
	return M_res;
}

Mat operator == (const Mat& M1, const Mat& M2)
{
	if(M1.data == NULL || M2.data == NULL)
	{
		cerr << "operator == (): Matrix has no data" << endl;
		exit(0);
	}
	if(M1.row != M2.row || M2.col != M2.col)
	{
		cerr << "operator == (): Matrix dimension do not match" << endl;
		exit(0);
	}
	Mat M_res(M1.row, M1.col);
	for(int i = 0; i < M1.row; i++)
	{
		for(int j = 0; j < M1.col; j++)
		{
			M_res.data[i][j] = (M1.data[i][j] == M2.data[i][j]);
		}
	}
	return M_res;
}

Mat::Mat(const Mat& M)
{
	if(M.data == NULL)
	{
		cerr << "Copy Mat(): Matrix has no data" << endl;
		exit(0);
	}
	row = M.row;
	col = M.col;
	data = new double*[row];
	for(int i = 0; i < row; i++)
	{
		data[i] = new double[col];
		for(int j = 0; j < col; j++)
		{
			data[i][j] = M.data[i][j];
		}
	}
}

Mat Mat::Diag()
{
	if(data == NULL)
	{
		cerr << "Diag(): Matrix has no data" << endl;
		exit(0);
	}
	int n = row > col ? col : row;
	Mat M_res(n, 1);
	for(int i = 0; i < n; i++)
	{
		M_res.data[i][0] = data[i][i];
	}
	return M_res;
}

Mat Mat::SVD(Mat& U, Mat& W, Mat& V)
{
	if(data == NULL)
	{
		cerr << "SVD(): Matrix has no data" << endl;
		exit(0);
	}
	if(U.row != row || U.col != row || W.row != row || W.col != col || V.row != col || V.col != col)
	{
		cerr << "SVD(): Matrix dimension do not match" << endl;
		exit(0);
	}
	double* a = new double[row * col];
	double* u = new double[row * row];
	double* v = new double[col * col];
	double* w = new double[row * col];

	for(int i = 0; i < row; i++)
	{
		for(int j = 0; j < col; j++)
		{
			a[i * col + j] = data[i][j];
		}
	}

	int ka = row > col ? row + 1 : col + 1;
	bmuav(a, row, col, u, v, 0.00001, ka);

	Mat M_res(row, col);

	for(i = 0; i < row; i++)
	{
		for(int j = 0; j < col; j++)
		{
			M_res.data[i][j] = a[i * col + j];
			W.data[i][j] = w[i * col + j];
		}
		for(j = 0; j < row; j++)
		{
			U.data[i][j] = u[i * row + j];
		}
	}

	for(i = 0; i < col; i++)
	{
		for(int j = 0; j < col; j++)
		{
			V.data[i][j] = v[i * col + j];
		}
	}

	delete[] w;
	delete[] v;
	delete[] u;
	delete[] a;
	return M_res;
}

Mat& Mat::operator += (const Mat& M)
{
	if(data == NULL || M.data == NULL)
	{
		cerr << "operator += (): Matrix has no data" << endl;
		exit(0);
	}
	if(row != M.row || col != M.col)
	{
		cerr << "operator += (): Matrix dimension do not match" << endl;
		exit(0);
	}
	for(int i = 0; i < row; i++)
	{
		for(int j = 0; j < col; j++)
		{
			data[i][j] += M.data[i][j];
		}
	}
	return *this;
}

Mat& Mat::operator -= (const Mat& M)
{
	if(data == NULL || M.data == NULL)
	{
		cerr << "operator -= (): Matrix has no data" << endl;
		exit(0);
	}
	if(row != M.row || col != M.col)
	{
		cerr << "operator -= (): Matrix dimension do not match" << endl;
		exit(0);
	}
	for(int i = 0; i < row; i++)
	{
		for(int j = 0; j < col; j++)
		{
			data[i][j] -= M.data[i][j];
		}
	}
	return *this;
}

Mat& Mat::operator *= (const double& s)
{
	if(data == NULL)
	{
		cerr << "operator *= (): Matrix has no data" << endl;
		exit(0);
	}
	for(int i = 0; i < row; i++)
	{
		for(int j = 0; j < col; j++)
		{
			data[i][j] *= s;
		}
	}
	return *this;
}

Mat& Mat::operator /= (const Mat& M)
{
	if(data == NULL || M.data == NULL)
	{
		cerr << "operator /= (): Matrix has no data" << endl;
		exit(0);
	}

	if(row != M.row || col != M.col)
	{
		cerr << "operator /= (): Matrix dimension do not match" << endl;
		exit(0);
	}

	for(int i = 0; i < row; i++)
	{
		for(int j = 0; j < col; j++)
		{
			if(M.data[i][j] == 0)
			{
				cerr << "operator /= (): Matrix has zero element" << endl;
				exit(0);
			}
			data[i][j] /= M.data[i][j];
		}
	}
	return *this;
}

Mat Mat::MSub(const int& li, const int& lj, const int& srow, const int& scol)
{
	if(data == NULL)
	{
		cerr << "MSub(): Matrix has no data" << endl;
		exit(0);
	}
	if(li < 0 || lj < 0 || li >= row || lj >= col)
	{
		cerr << "MSub(): index exceed matrix dimension" << endl;
		exit(0);
	}
	if(srow <= 0 || scol <=0 || li + srow > row || lj + scol > col)
	{
		cerr << "MSub(): sub size exceed matrix size" << endl;
		exit(0);
	}
	/*Mat* pM_res = new Mat(srow, scol);
	for(int i = 0; i < srow; i++)
	{
	for(int j = 0; j < scol; j++)
	{
	pM_res->data[i][j] = data[i + li][j + lj];
	}
	}
	return *pM_res;*/
	Mat M_res(srow, scol);
	for(int i = 0; i < srow; i++)
	{
		for(int j = 0; j < scol; j++)
		{
			M_res.data[i][j] = data[i + li][j + lj];
		}
	}
	return M_res;
}

void Mat::MPub(const Mat& M, const int& li, const int& lj)
{
	if(data == NULL || M.data == NULL)
	{
		cerr << "MPub(): Matrix has no data" << endl;
		exit(0);
	}
	if(li < 0 || lj < 0 || li >= row || lj >= col)
	{
		cerr << "MPub(): index exceed matrix dimension" << endl;
		exit(0);
	}
	if(li + M.row > row || lj + M.col > col)
	{
		cerr << "MPub(): pub size exceed matrix size" << endl;
		exit(0);
	}
	for(int i = 0; i < M.row; i++)
	{
		for(int j = 0; j < M.col; j++)
		{
			data[i + li][j + lj] = M.data[i][j];
		}
	}
}

Mat Mat::operator - ()
{
	if(data == NULL)
	{
		cerr << "operator -(): Matrix has no data" << endl;
		exit(0);
	}
	Mat M_res(row, col);
	for(int i = 0; i < row; i++)
	{
		for(int j = 0; j < col; j++)
		{
			M_res.data[i][j] = -data[i][j];
		}
	}
	return M_res;
}

Mat MRepmat(const Mat& M, const int& m, const int& n)
{
	if(M.data == NULL)
	{
		cerr << "MRepmat(): Matrix has no data" << endl;
		exit(0);
	}
	if(n <= 0)
	{
		cerr << "MRepmat(): n <= 0" << endl;
		exit(0);
	}
	int n_row = M.row * m;
	int n_col = M.col * n;
	
	Mat M_res(n_row, n_col);
	for(int i = 0; i < m; i++)
	{
		for(int j = 0; j < n; j++)
		{
			for(int p = 0; p < M.row; p++)
			{
				for(int q = 0; q < M.col; q++)
				{
					M_res.data[i * M.row + p][j * M.col + q] = M.data[p][q];
				}
			}
		}
	}
	return M_res;
}

Mat MMean(const Mat& M)
{
	if(M.data == NULL)
	{
		cerr << "MMean(): Matrix has no data" << endl;
		exit(0);
	}
	Mat M_res(1, M.col);
	for(int i = 0; i < M_res.col; i++)
	{
		double sum = 0;
		for(int j = 0; j < M.row; j++)
		{
			sum += M.data[j][i];
		}
		sum /= M.row;
		M_res.data[0][i] = sum;
	}
	return M_res;
}

Mat MDPow(const Mat& M, const double& q)
{
	if(M.data == NULL)
	{
		cerr << "MDPow(): Matrix has no data" << endl;
		exit(0);
	}
	Mat M_res(M.row, M.col);
	for(int i = 0; i < M.row; i++)
	{
		for(int j = 0; j < M.col; j++)
		{
			M_res.data[i][j] = pow(M.data[i][j], q);
		}
	}
	return M_res;
}

Mat MSum(const Mat& M, int type = 1)
{
	if(M.data == NULL)
	{
		cerr << "MSum(): Matrix has no data" << endl;
		exit(0);
	}
	if(type != 1 && type != 2)
	{
		cerr << "MSum(): type error" << endl;
		exit(0);
	}
	if(type == 1)
	{
		Mat M_res(1, M.col);
		for(int i = 0; i < M.col; i++)
		{
			double sum = 0;
			for(int j = 0; j < M.row; j++)
			{
				sum += M.data[j][i];
			}
			M_res.data[0][i] = sum;
		}
		return M_res;
	}
	else
	{
		Mat M_res(M.row, 1);
		for(int i = 0; i < M.row; i++)
		{
			double sum = 0;
			for(int j = 0; j < M.col; j++)
			{
				sum += M.data[i][j];
			}
			M_res.data[i][0] = sum;
		}
		return M_res;
	}
}

double** MReadFile(const char filename[], const int& row, const int& col)
{
	if(row <= 0 || col <= 0)
	{
		cerr << "MReadFile(): row or col <= 0" << endl;
		exit(0);
	}
	if(filename == 0)
	{
		cerr << "MReadFile(): filename error" << endl;
		exit(0);
	}
	ifstream in(filename);
	cout << "Reading file ... " << endl;
	double** data = new double*[row];
	for(int i = 0; i < row; i++)
	{
		cout << (i + 1) * 100 / row << " %\r" << flush;
		data[i] = new double[col];
		for(int j = 0; j < col; j++)
		{
			in >> data[i][j];
		}
	}
	in.close();
	cout << "Reading file finish." << endl;
	return data;
}

double** MReadFileBin(const char filename[], const int& row, const int& col)
{
	if(row <= 0 || col <= 0)
	{
		cerr << "MReadFile(): row or col <= 0" << endl;
		exit(0);
	}
	if(filename == 0)
	{
		cerr << "MReadFile(): filename error" << endl;
		exit(0);
	}
	ifstream in(filename, ios::binary);
	cout << "Reading file ... " << endl;
	double** data = new double*[row];
	for(int i = 0; i < row; i++)
	{
		cout << (i + 1) * 100 / row << " %\r" << flush;
		data[i] = new double[col];
		for(int j = 0; j < col; j++)
		{
			in >> data[i][j];
		}
	}
	in.close();
	cout << "Reading file finish." << endl;
	return data;
}

void MWriteFile(const char filename[], const Mat& M)
{
	if(filename == 0)
	{
		cerr << "MWriteFile(): filename error" << endl;
		exit(0);
	}
	ofstream out(filename);
	cout << "Writing file ... " << endl;
	for(int i = 0; i < M.row; i++)
	{
		cout << (i + 1) * 100 / M.row << " %\r" << flush;
		for(int j = 0; j < M.col; j++)
		{
			out << M.data[i][j] << "\t";
		}
		out << endl;
	}
	out.close();
	cout << "Writing file finish." << endl;
}

void MWriteFileBin(const char filename[], const Mat& M)
{
	if(filename == 0)
	{
		cerr << "MWriteFile(): filename error" << endl;
		exit(0);
	}
	ofstream out(filename, ios::binary);
	cout << "Writing file ... " << endl;
	for(int i = 0; i < M.row; i++)
	{
		cout << (i + 1) * 100 / M.row << " %\r" << flush;
		for(int j = 0; j < M.col; j++)
		{
			out << M.data[i][j] << "\t";
		}
		out << endl;
	}
	out.close();
	cout << "Writing file finish." << endl;
}

Mat::Mat(const int& row, const int& col, double** d_data)
{
	if(row <= 0 || col <= 0)
	{
		cerr << "Mat(): row or col <= 0" << endl;
		exit(0);
	}
	if(d_data == NULL)
	{
		cerr << "Mat(): Data is NULL" << endl;
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
			data[i][j] = d_data[i][j];
		}
	}
}
Mat MEye(const int& n)
{
	if(n <= 0)
	{
		cerr << "MEye(): n <= 0" << endl;
		exit(0);
	}
	Mat M_res(n, n);
	for(int i = 0; i < n; i++)
	{
		M_res.data[i][i] = 1;
	}
	return M_res;
}

double* & Mat::operator [] (const int& index)
{
	if(data == NULL)
	{
		cerr << "operator [] (): Matrix has no data" << endl;
		exit(0);
	}
	if(index < 0 || index >= row)
	{
		cerr << "operator [] (): Index exceed dimension" << endl;
		exit(0);
	}
	return data[index];
}
Mat::Mat(const int& row, const int& col)
{
	if(row <= 0 || col <= 0)
	{
		cerr << "Mat(): row or col <= 0" << endl;
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
	//cout << "Unconstruction !" << endl;
	if(data == NULL)
	{
		cerr << "~Mat(): data == NULL" << endl;
		exit(0);
	}
	for(int i = 0; i < row; i++)
	{
		if(data[i] == NULL)
		{
			cerr << "~Mat(): data == NULL" << endl;
			exit(0);
		}
		delete[] data[i];
	}
	delete[] data;
}
void Mat::MShow()
{
	if(data == NULL)
	{
		cerr << "MShow(): Matrix has no data" << endl;
		exit(0);
	}
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
		cerr << "operator = (): Matrix has no data" << endl;
		exit(0);
	}
	if(M.row != row || M.col != col)
	{
		cerr << "operator = (): Matrix dimension do not match" << endl;
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
Mat MInv(const Mat& M)
{
	if(M.data == NULL)
	{
		cerr << "MInv(): Matrix has no data" << endl;
		exit(0);
	}
	if(M.col != M.row)
	{
		cerr << "MInv(): Matrix is not square" << endl;
		exit(0);
	}
	double d_value;
	int i_n = M.row;
	
	Mat M_res(i_n, i_n);
	M_res = M;
	
	for(int i = 0; i < i_n; i++)
	{
		d_value = 1.0 / M_res.data[i][i];
		M_res.data[i][i] = d_value;
		for(int j = 0; j < i_n; j++)
		{
			if(i != j)
			{
				M_res.data[i][j] *= -d_value;
			}
		}
		for(j = 0; j < i_n; j++)
		{
			if(i != j)
			{
				M_res.data[j][i] *= d_value;
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
						M_res.data[j][k] += M_res.data[j][i] * M_res.data[i][k] / d_value;
					}
				}
			}
		}
	}
	return M_res;
}

Mat operator +(const Mat& M1, const Mat& M2)
{
	if(M1.data == NULL || M2.data == NULL)
	{
		cerr << "operator + (): M1 or M2 has no data" << endl;
		exit(0);
	}
	if(M1.row != M2.row || M1.col != M2.col)
	{
		cerr << "operator + (): Matrix dimension do not match" << endl;
		exit(0);
	}
	Mat M_res(M1.row, M1.col);
	for(int i = 0; i < M_res.row; i++)
	{
		for(int j = 0; j < M_res.col; j++)
		{
			M_res.data[i][j] = M1.data[i][j] + M2.data[i][j];
		}
	}
	return M_res;
}

Mat operator -(const Mat& M1, const Mat& M2)
{
	if(M1.data == NULL || M2.data == NULL)
	{
		cerr << "operator - (): M1 or M2 has no data" << endl;
		exit(0);
	}
	if(M1.row != M2.row || M1.col != M2.col)
	{
		cerr << "operator - (): Matrix dimension do not match" << endl;
		exit(0);
	}
	Mat M_res(M1.row, M1.col);
	for(int i = 0; i < M_res.row; i++)
	{
		for(int j = 0; j < M_res.col; j++)
		{
			M_res.data[i][j] = M1.data[i][j] - M2.data[i][j];
		}
	}
	return M_res;
}

Mat operator -(const double& s, const Mat& M2)
{
	if(M2.data == NULL)
	{
		cerr << "operator - (): M2 has no data" << endl;
		exit(0);
	}
	Mat M_res(M2.row, M2.col);
	for(int i = 0; i < M_res.row; i++)
	{
		for(int j = 0; j < M_res.col; j++)
		{
			M_res.data[i][j] = s - M2.data[i][j];
		}
	}
	return M_res;
}

Mat operator -(const Mat& M, const double& s)
{
	if(M.data == NULL)
	{
		cerr << "operator - (M, s): M has no data" << endl;
		exit(0);
	}
	Mat M_res(M.row, M.col);
	for(int i = 0; i < M_res.row; i++)
	{
		for(int j = 0; j < M_res.col; j++)
		{
			M_res.data[i][j] = M.data[i][j] - s;
		}
	}
	return M_res;
}

Mat operator *(const Mat& M1, const Mat& M2)
{
	if(M1.data == NULL || M2.data == NULL)
	{
		cerr << "operator * (): M1 or M2 has no data" << endl;
		exit(0);
	}
	if(M1.col != M2.row)
	{
		cerr << "operator * (): Matrix dimension do not match" << endl;
		exit(0);
	}
	/*Mat* pM_res = new Mat(M1.row, M2.col);
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
	return *pM_res;*/
	Mat M_res(M1.row, M2.col);
	for(int i = 0; i < M_res.row; i++)
	{
		for(int j = 0; j < M_res.col; j++)
		{
			for(int k = 0; k < M1.col; k++)
			{
				M_res.data[i][j] += M1.data[i][k] * M2.data[k][j];
			}
		}
	}
	return M_res;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
Mat operator *(const double& s, const Mat& M2)
{
	if(M2.data == NULL)
	{
		cerr << "operator * (): M2 has no data" << endl;
		exit(0);
	}
	
	/*
	Mat* pM_res = new Mat(M2.row, M2.col);
	for(int i = 0; i < pM_res->row; i++)
	{
	for(int j = 0; j < pM_res->col; j++)
	{
	pM_res->data[i][j] = s * M2.data[i][j];
	}
	}
	return *pM_res;*/
	
	Mat M_res(M2.row, M2.col);
	for(int i = 0; i < M_res.row; i++)
	{
		for(int j = 0; j < M_res.col; j++)
		{
			M_res.data[i][j] = s * M2.data[i][j];
		}
	}
	return M_res;
}

Mat MTrans(const Mat& M)
{
	if(M.data == NULL)
	{
		cerr << "MTrans(): Matrix has no data" << endl;
		exit(0);
	}
	Mat M_res(M.col, M.row);
	for(int i = 0; i < M_res.row; i++)
	{
		for(int j = 0; j < M_res.col; j++)
		{
			M_res.data[i][j] = M.data[j][i];
		}
	}
	return M_res;
}

Mat MSigmoid(const Mat& M)
{
	if(M.data == NULL)
	{
		cerr << "MSigmoid(): Matrix has no data" << endl;
		exit(0);
	}
	Mat M_res(M.row, M.col);
	for(int i = 0; i < M_res.row; i++)
	{
		for(int j = 0; j < M_res.col; j++)
		{
			M_res.data[i][j] = 1.0 / (1.0 + exp(-M.data[i][j]));
		}
	}
	return M_res;
}

Mat MDiag(const Mat& M)
{
	if(M.data == NULL)
	{
		cerr << "MDiag(): Matrix has no data" << endl;
		exit(0);
	}
	if(M.col != 1)
	{
		cerr << "MDiag(): Matrix is not cvector" << endl;
		exit(0);
	}
	Mat Mat_res(M.row, M.row);
	for(int i = 0; i < Mat_res.row; i++)
	{
		Mat_res.data[i][i] = M.data[i][0];
	}
	return Mat_res;
}

bool Mat::QR(Mat& Q, Mat& R)
{
	if(data == NULL)
	{
		cerr << "QR(): Matrix has no data" << endl;
		exit(0);
	}
	if(Q.data == NULL || R.data == NULL)
	{
		cerr << "QR(): Q or R has not init data" << endl;
		exit(0);
	}
	if(row != col)
	{
		cerr << "QR(): Matrix is not square" << endl;
		exit(0);
	}
	if(Q.row != Q.col || R.row != R.col || Q.row != row)
	{
		cerr << "QR(): Matrix dimension do not match" << endl;
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

Mat MReshape(const Mat& M, const int& row, const int& col)
{
	if(M.data == NULL)
	{
		cerr << "MReshape(): Matrix has no data" << endl;
		exit(0);
	}
	if(row <=0 || col <= 0)
	{
		cerr << "MReshape(): row or col <= 0" << endl;
		exit(0);
	}
	if(M.row * M.col != row * col)
	{
		cerr << "MReshape(): Matrix dimension do not match" << endl;
		exit(0);
	}
	Mat M_res(row, col);
	for(int i = 0; i < col; i++)
	{
		for(int j = 0; j < row; j++)
		{
			int index = i * row + j;
			int oi = index / M.row;
			int oj = index % M.row;
			M_res.data[j][i] = M.data[oj][oi];
		}
	}
	return M_res;
}

Mat operator / (const Mat& M, const double& s)
{
	if(M.data == NULL)
	{
		cerr << "operator / (): Matrix has no data" << endl;
		exit(0);
	}
	if(s == 0)
	{
		cerr << "operator / (): s = 0" << endl;
		exit(0);
	}
	Mat M_res(M.row, M.col);
	for(int i = 0; i < M.row; i++)
	{
		for(int j = 0; j < M.col; j++)
		{
			M_res.data[i][j] = M.data[i][j] / s;
		}
	}
	return M_res;
}

Mat operator / (const Mat& M, const Mat& s)
{
	if(M.data == NULL || s.data == NULL)
	{
		cerr << "operator / (): Matrix has no data" << endl;
		exit(0);
	}
	if(s.row != 1 || s.col != 1)
	{
		cerr << "operator / (): s is not a scalar" << endl;
		exit(0);

	}
	if(s.data[0][0] == 0)
	{
		cerr << "operator / (): s = 0" << endl;
		exit(0);
	}
	Mat M_res(M.row, M.col);
	for(int i = 0; i < M.row; i++)
	{
		for(int j = 0; j < M.col; j++)
		{
			M_res.data[i][j] = M.data[i][j] / s.data[0][0];
		}
	}
	return M_res;
}

Mat operator / (const double& s, const Mat& M)
{
	if(M.data == NULL)
	{
		cerr << "operator / (): Matrix has no data" << endl;
		exit(0);
	}
	Mat M_res(M.row, M.col);
	for(int i = 0; i < M.row; i++)
	{
		for(int j = 0; j < M.col; j++)
		{
			if(M.data[i][j] == 0)
			{
				cerr << "operator / (): Matrix has zero element" << endl;
				exit(0);
			}
			M_res.data[i][j] = s / M.data[i][j];
		}
	}
	return M_res;
}

Mat operator * (const Mat& M, double s)
{
	if(M.data == NULL)
	{
		cerr << "operator * (): Matrix has no data" << endl;
		exit(0);
	}
	Mat M_res(M.row, M.col);
	for(int i = 0; i < M.row; i++)
	{
		for(int j = 0; j < M.col; j++)
		{
			M_res.data[i][j] = M.data[i][j] * s;
		}
	}
	return M_res;
}

Mat operator + (const Mat& M, double s)
{
	if(M.data == NULL)
	{
		cerr << "operator + (): Matrix has no data" << endl;
		exit(0);
	}
	Mat M_res(M.row, M.col);
	for(int i = 0; i < M.row; i++)
	{
		for(int j = 0; j < M.col; j++)
		{
			M_res.data[i][j] = M.data[i][j] + s;
		}
	}
	return M_res;
}

Mat MMax(const Mat& M, double s)
{
	if(M.data == NULL)
	{
		cerr << "MMax(): Matrix has no data" << endl;
		exit(0);
	}
	Mat M_res(M.row, M.col);
	for(int i = 0; i < M.row; i++)
	{
		for(int j = 0; j < M.col; j++)
		{
			if(M.data[i][j] > s)
			{
				M_res.data[i][j] = M.data[i][j];
			}
			else
			{
				M_res.data[i][j] = s;
			}
		}
	}
	return M_res;
}

Mat MMin(const Mat& M, double s)
{
	if(M.data == NULL)
	{
		cerr << "MMin(): Matrix has no data" << endl;
		exit(0);
	}
	Mat M_res(M.row, M.col);
	for(int i = 0; i < M.row; i++)
	{
		for(int j = 0; j < M.col; j++)
		{
			if(M.data[i][j] < s)
			{
				M_res.data[i][j] = M.data[i][j];
			}
			else
			{
				M_res.data[i][j] = s;
			}
		}
	}
	return M_res;
}

Mat MStd(const Mat& M)
{
	if(M.data == NULL)
	{
		cerr << "MStd(): Matrix has no data" << endl;
		exit(0);
	}
	Mat M_res(1, M.col);
	M_res = MDPow(MSum(MDPow(M - MRepmat(MMean(M), M.row, 1), 2)) / (M.row - 1), 0.5);
	
	return M_res;
}

Mat MLog(const Mat& M)
{
	if(M.data == NULL)
	{
		cerr << "MLog(): Matrix has no data" << endl;
		exit(0);
	}
	Mat M_res(M.row, M.col);
	for(int i = 0; i < M.row; i++)
	{
		for(int j = 0; j < M.col; j++)
		{
			if(M.data[i][j] <= 0)
			{
				cerr << "MLog(): Matrix has negative element" << endl;
				exit(0);
			}
			M_res.data[i][j] = log(M.data[i][j]);
		}
	}
	return M_res;
}

Mat MDProduct(const Mat& M1, const Mat& M2)
{
	if(M1.data == NULL || M2.data == NULL)
	{
		cerr << "MDProduct(): Matrix has no data" << endl;
		exit(0);
	}
	if(M1.row != M2.row || M1.col != M2.col)
	{
		cerr << "MDProduct(): Matrix dimension do not match" << endl;
		exit(0);
	}
	
	Mat M_res(M1.row, M1.col);
	
	for(int i = 0; i < M1.row; i++)
	{
		for(int j = 0; j < M1.col; j++)
		{
			M_res.data[i][j] = M1.data[i][j] * M2.data[i][j];
		}
	}
	
	return M_res;
}

Mat MSigmoidInv(const Mat& M)
{
	if(M.data == NULL)
	{
		cerr << "MSigmoidInv(): Matrix has no data" << endl;
		exit(0);
	}
	
	Mat M_res(M.row, M.col);
	M_res = MDProduct(MSigmoid(M), (1 - MSigmoid(M)));
	
	return M_res;
}

Mat MOnes(const int& m, const int& n)
{
	if(m == 0 || n == 0)
	{
		cerr << "MOnes(): m or n equals zero" << endl;
		exit(0);
	}
	Mat M_res(m, n);
	for(int i = 0; i < m; i++)
	{
		for(int j = 0; j < n; j++)
		{
			M_res.data[i][j] = 1;
		}
	}
	return M_res;
}

Mat MDDivide(const Mat& M1, const Mat& M2)
{
	if(M1.data == NULL || M2.data == NULL)
	{
		cerr << "MDDivide(): Matrix M1 or M2 has no data" << endl;
		exit(0);
	}
	if(M1.row != M2.row || M1.col != M2.col)
	{
		cerr << "MDDivide(): Matrix dimension do not match" << endl;
		exit(0);
	}
	Mat M_res(M1.row, M1.col);
	for(int i = 0; i < M_res.row; i++)
	{
		for(int j = 0; j < M_res.col; j++)
		{
			if(M2.data[i][j] == 0)
			{
				cerr << "MDDivide(): M2 has zero element" << endl;
				exit(0);
			}
			M_res.data[i][j] = M1.data[i][j] / M2.data[i][j];
		}
	}
	return M_res;
}

Mat operator <= (const Mat& M, const double& s)
{
	if(M.data == NULL)
	{
		cerr << "operator <= (): Matrix has no data" << endl;
		exit(0);
	}
	Mat M_res(M.row, M.col);
	for(int i = 0; i < M.row; i++)
	{
		for(int j = 0; j < M.col; j++)
		{
			if(M.data[i][j] <= s)
			{
				M_res.data[i][j] = 1;
			}
		}
	}
	return M_res;
}

Mat MCusum(const Mat& M)
{
	if(M.data == NULL)
	{
		cerr << "MCusum(): Matrix has no data" << endl;
		exit(0);
	}
	Mat M_res(M.row, M.col);
	for(int i = 0; i < M.row; i++)
	{
		for(int j = 0; j < M.col; j++)
		{
			if(i == 0)
			{
				M_res.data[i][j] = M.data[i][j];
			}
			else
			{
				M_res.data[i][j] = M_res.data[i - 1][j] + M.data[i][j];
			}
		}
	}
	return M_res;
}

Mat MAbs(const Mat& M)
{
	if(M.data == NULL)
	{
		cerr << "MAbs(): Matrix has no data" << endl;
		exit(0);
	}
	Mat M_res(M.row, M.col);
	for(int i = 0; i < M.row; i++)
	{
		for(int j = 0; j < M.col; j++)
		{
			M_res.data[i][j] = fabs(M.data[i][j]);
		}
	}
	return M_res;
}

Mat Max(const Mat& M, Mat& order, int type = 1)
{
	if(M.data == NULL)
	{
		cerr << "Max(): Matrix has no data" << endl;
		exit(0);
	}
	Mat M_res(1, M.col);
	for(int i = 0; i < M.col; i++)
	{
		double dm = M.data[0][i];
		for(int j = 0; j < M.row; j++)
		{
			if(dm < M.data[j][i])
			{
				dm = M.data[j][i];
				order[0][i] = j;
			}
		}
		M_res.data[0][i] = dm;
	}
	return M_res;
}

Mat Max(const Mat& M, int type = 1)
{
	if(M.data == NULL)
	{
		cerr << "Max(): Matrix has no data" << endl;
		exit(0);
	}
	Mat M_res(1, M.col);
	for(int i = 0; i < M.col; i++)
	{
		double dm = M.data[0][i];
		for(int j = 0; j < M.row; j++)
		{
			if(dm < M.data[j][i])
			{
				dm = M.data[j][i];
			}
		}
		M_res.data[0][i] = dm;
	}
	return M_res;
}

Mat MExp(const Mat& M)
{
	if(M.data == NULL)
	{
		cerr << "MExp(): Matrix has no data" << endl;
		exit(0);
	}
	Mat M_res(M.row, M.col);
	for(int i = 0; i < M.row; i++)
	{
		for(int j = 0; j < M.col; j++)
		{
			M_res.data[i][j] = exp(M.data[i][j]);
			if(M_res.data[i][j] <= 0)
			{
				cout << "!" << endl;
				cout << M.data[i][j] << endl;
				system("pause");
			}
		}
	}
	return M_res;
}


#endif
