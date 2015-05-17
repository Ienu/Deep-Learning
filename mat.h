// MAT_H
#ifndef MAT_H
#define MAT_H

#include <iostream>

using namespace std;

class Mat
{
public:
	double** data;
	
public:
	int row;
	int col;
	
public:
	
	Mat& operator =(const Mat& M);	
	void MShow();	
	Mat(int row, int col);
	~Mat();
	
};

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
Mat MProduct(const Mat& M1, const Mat& M2)
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
#endif
