//minFunc_H
#ifndef MINFUNC_H
#define MINFUNC_H

#include "mat.h"

#define LIMITL	10

// type 1 BFGS
// type 2 L-BFGS

struct minFuncOp
{
	int type;
	int maxIter;
};

Mat minFunc(
			double (*fun)
			(Mat&, 
			const int&, 
			const int&, 
			const double&, 
			const double&, 
			const double&, 
			const Mat&, 
			Mat&),
			Mat& theta, 
			const int& visibleSize, 
			const int& hiddenSize,
			const double& lambda,
			const double& sparsityParam,
			const double& beta,
			const Mat& data,
			Mat& grad,
			const minFuncOp& option)
{
	int n = 2 * visibleSize * hiddenSize + visibleSize + hiddenSize;
	
	Mat g(n, 1);
	fun(theta, visibleSize, hiddenSize, lambda, sparsityParam, beta, data, g);
	
	Mat old_dirs(n, option.maxIter);
	Mat old_stps(n, option.maxIter);
	
	//Mat H(n, n);
	Mat d(n, 1);
	Mat y(n, 1);
	Mat g_old(n, 1);
	Mat s(n, 1);
	Mat ys(1, 1);
	
	Mat x(n, 1);
	x = theta;
	
	time_t t_start = time(NULL);
	time_t t_end;
	
	//double cost;
	double t = 1.0;
	int numCorrections = 0;
	
	for(int i = 0; i < option.maxIter; i++)
	{
		t_end = time(NULL);
		cout << "Iter : " << i << " last for " << difftime(t_end, t_start) << " s" << endl;
		switch(option.type)
		{
		/*case 0:	// gradient descend
		M_xNew = M_xOld - M_gOld;
		fun(M_xNew, visibleSize, hiddenSize, lambda, sparsityParam, beta, data, M_gNew);
		
		  M_gOld = M_gNew;
		  M_xOld = M_xNew;
			break;*/
		case 2:	// L-BFGS
			if(i == 0)
			{
				d = -g;
				//H = MEye(n);
			}
			else
			{
				y = g - g_old;
				s = t * d;
				ys = MTrans(y) * s;
				
				old_dirs.MPub(s, 0, numCorrections);
				old_stps.MPub(y, 0, numCorrections);
				numCorrections++;

				int begin = 0;
				int k = 0;

				if(numCorrections >= LIMITL)
				{
					begin = numCorrections - LIMITL;
					k = LIMITL;
				}
				else
				{
					k = numCorrections;
					//cerr << "Not deal yet" << endl;
					//exit(0);
				}
				double h = ys.data[0][0] / (MTrans(y) * y).data[0][0];
				
				

				Mat rou(1, k);
				rou = MDDivide(MOnes(1, k), MSum(MDProduct(old_stps.MSub(0, begin, n, k), 
					old_dirs.MSub(0, begin, n, k))));// + 1e-9 * MOnes(1, k));

				Mat q(n, k + 1);
				Mat r(n, k + 1);
				Mat Alpha(k, 1);
				Mat Beta(k, 1);
				
				q.MPub(-g, 0, k);
				for(int j = k - 1; j >= 0; j--)
				{
					Alpha[j][0] = rou[0][j] * (MTrans(old_dirs.MSub(0, j + begin, n, 1)) * q.MSub(0, j + 1, n, 1)).data[0][0];
					q.MPub(q.MSub(0, j + 1, n, 1) - Alpha[j][0] * old_stps.MSub(0, j + begin, n, 1), 0, j);
				}
				r.MPub(h * q.MSub(0, 0, n, 1), 0, 0);

				for(j = 0; j < k; j++)
				{
					Beta[j][0] = rou[0][j] * (MTrans(old_stps.MSub(0, j + begin, n, 1)) * r.MSub(0, j, n, 1)).data[0][0];
					r.MPub(r.MSub(0, j, n, 1) + old_dirs.MSub(0, j + begin, n, 1) * (Alpha[j][0] - Beta[j][0]), 0, j + 1);
				}
				d = r.MSub(0, k, n, 1);
			}
			g_old = g;
			if(i == 0)
			{
				double dm = 1.0 / (MSum(MAbs(g)).data[0][0]);// + 1e-9);
				t = 1.0 > dm ? dm : 1.0;
			}
			else
			{
				t = 1;
			}
			x = x + t * d;
			fun(x, visibleSize, hiddenSize, lambda, sparsityParam, beta, data, g);
			break;
		}
	}
	return x;
}



#endif
