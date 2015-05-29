#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include "mat.h"

#include "minfunc.h"

using namespace std;

const int patchSize = 8;
const int numPatches = 10000;

const int visibleSize = patchSize * patchSize;
const int hiddenSize = 25;

const double sparsityParam = 0.01;
const double lambda = 0.0001;

const double beta = 3;

Mat sampleIMAGES()
{
	int imageCount = 10;

	int row = 512;
	int col = 512;

	double** IMAGES = MReadFile("IMAGES.dat", imageCount * row, col);
	Mat M_patches(patchSize * patchSize, numPatches);

	srand(time(NULL));
	for(int i = 0; i < numPatches; i++)
	{
		int imageNum = rand() % imageCount;
		int xPos = rand() % (col - patchSize);
		int yPos = rand() % (row - patchSize);

		for(int m = 0; m < patchSize; m++)
		{
			for(int n = 0; n < patchSize; n++)
			{
				M_patches.data[m * patchSize + n][i] 
					= IMAGES[imageNum * row + yPos + n][xPos + m];
			}
		}
	}
	M_patches = M_patches -(MSum(MSum(M_patches), 2).data[0][0]) / double(patchSize * patchSize);

	Mat M_pstd(1, 1);
	M_pstd = MStd(MReshape(M_patches, patchSize * patchSize * numPatches, 1));
	
	double pstd = 3 * M_pstd[0][0];
	cout << pstd << endl;

	M_patches = MMax(MMin(M_patches, pstd), -pstd) / pstd * 0.4 + 0.5;

	return M_patches;
}

Mat initializeParameters(int hiddenSize, int visibleSize)
{
	double r = sqrt(6.0 / (hiddenSize + visibleSize + 1));
	Mat M_theta(hiddenSize * visibleSize * 2 + hiddenSize + visibleSize, 1);

	srand(time(NULL));
	for(int i = 0; i < hiddenSize * visibleSize * 2; i++)
	{
		M_theta.data[i][0] = (rand() % 10000) / 10000.0 * 2 * r - r;
	}

	return M_theta;
}

double sparseAutoEncoderCost
(Mat& theta, 
 const int& visibleSize, 
 const int& hiddenSize,
 const double& lambda,
 const double& sparsityParam,
 const double& beta,
 const Mat& data,
 Mat& grad)
{
	Mat W1(hiddenSize, visibleSize);
	Mat W2(visibleSize, hiddenSize);

	Mat b1(hiddenSize, 1);
	Mat b2(visibleSize, 1);

	W1 = MReshape(theta.MSub(0, 0, hiddenSize * visibleSize, 1), hiddenSize, visibleSize);
	W2 = MReshape(theta.MSub(hiddenSize * visibleSize, 0, hiddenSize * visibleSize, 1), visibleSize, hiddenSize);
	b1 = theta.MSub(hiddenSize * visibleSize * 2, 0, hiddenSize, 1);
	b2 = theta.MSub(hiddenSize * visibleSize * 2 + hiddenSize, 0, visibleSize, 1);

	Mat cost(1, 1);

	Mat W1grad(hiddenSize, visibleSize);
	Mat W2grad(visibleSize, hiddenSize);

	Mat b1grad(hiddenSize, 1);
	Mat b2grad(visibleSize, 1);

	//Mat Jcost(1, 1);
	//Mat Jweight(1, 1);
	//Mat Jsparse(1, 1);

	Mat z2(hiddenSize, numPatches);
	z2 = W1 * data + MRepmat(b1, 1, numPatches);

	//MWriteFile("z2c.dat", z2);

	Mat a2(hiddenSize, numPatches);
	a2 = MSigmoid(z2);

	Mat z3(visibleSize, numPatches);
	z3 = W2 * a2 + MRepmat(b2, 1, numPatches);

	Mat a3(visibleSize, numPatches);
	a3 = MSigmoid(z3);

	//Jcost = (0.5 / numPatches) * MSum(MSum(MDPow(a3 - data, 2)), 2);
	//Jweight = 0.5 * (MSum(MSum(MDPow(W1, 2)), 2) + MSum(MSum(MDPow(W2, 2)), 2));

	Mat rho(hiddenSize, 1);
	rho = (1.0 / numPatches) * MSum(a2, 2);

	//Jsparse = MSum(sparsityParam * MLog(sparsityParam / rho) 
	//	+ (1 - sparsityParam) * MLog((1 - sparsityParam) / (1 - rho)));

	//cost = Jcost + lambda * Jweight + beta * Jsparse;

	Mat d3(visibleSize, numPatches);
	d3 = -MDProduct(data - a3, MSigmoidInv(z3));

	Mat sterm(hiddenSize, 1);
	sterm = beta * (-sparsityParam / rho + (1 - sparsityParam) / (1 - rho));

	Mat d2(hiddenSize, numPatches);
	d2 = MDProduct(MTrans(W2) * d3 + MRepmat(sterm, 1, numPatches), MSigmoidInv(z2));

	W1grad += d2 * MTrans(data);
	W1grad = (1.0 / numPatches) * W1grad + lambda * W1;

	W2grad += d3 * MTrans(a2);
	W2grad = (1.0 / numPatches) * W2grad + lambda * W2;

	b1grad += MSum(d2, 2);
	b1grad *= (1.0 / numPatches);

	b2grad += MSum(d3, 2);
	b2grad *= (1.0 / numPatches);

	grad.MPub(MReshape(W1grad, hiddenSize * visibleSize, 1), 0, 0);
	grad.MPub(MReshape(W2grad, hiddenSize * visibleSize, 1), hiddenSize * visibleSize, 0);
	grad.MPub(b1grad, hiddenSize * visibleSize * 2, 0);
	grad.MPub(b2grad, hiddenSize * visibleSize * 2 + hiddenSize, 0);

	return cost[0][0];
}

int main()
{

	time_t t_start;
	time_t t_end;
	
	double** patches = MReadFile("patches.dat", patchSize * patchSize, numPatches);

	Mat M_patches(patchSize * patchSize, numPatches, patches);

	//M_patches = sampleIMAGES();
	
	//double** W1 = MReadFile("W1.dat", hiddenSize, visibleSize);
	//double** W2 = MReadFile("W2.dat", visibleSize, hiddenSize);

	//MWriteFile("patches.dat", M_patches);
	
	double** th = MReadFile("THETA.dat", visibleSize * hiddenSize * 2 + visibleSize + hiddenSize, 1);

	Mat M_theta(hiddenSize * visibleSize * 2 + hiddenSize + visibleSize, 1, th);

	//M_theta = initializeParameters(hiddenSize, visibleSize);

	//MWriteFile("THETA.dat", M_theta);
	
	//Mat W1(hiddenSize, visibleSize);
	//Mat W2(visibleSize, hiddenSize);

	//W1 = MReshape(M_theta.MSub(0, 0, hiddenSize * visibleSize, 1), hiddenSize, visibleSize);
	//MWriteFile("W1.dat", W1);

	//W2 = MReshape(M_theta.MSub(hiddenSize * visibleSize, 0, hiddenSize * visibleSize, 1), visibleSize, hiddenSize);
	//MWriteFile("W2.dat", W2);
	
	//system("pause");

	Mat M_grad(hiddenSize * visibleSize * 2 + hiddenSize + visibleSize, 1);

	t_start = time(NULL);
	double cost = sparseAutoEncoderCost
		(M_theta, 
		visibleSize, 
		hiddenSize,
		lambda,
		sparsityParam,
		beta,
		M_patches,
		M_grad);

	minFuncOp opt;
	opt.maxIter = 400;
	opt.type = 2;

	Mat M_optTheta(hiddenSize * visibleSize * 2 + hiddenSize + visibleSize, 1);

	M_optTheta = minFunc(&sparseAutoEncoderCost, M_theta, visibleSize, hiddenSize, lambda, sparsityParam, beta, M_patches, M_grad, opt);
	MWriteFile("res.dat", M_optTheta);

	t_end = time(NULL);

	cout << "time last for " << difftime(t_end, t_start) << "s." << endl;

	Mat W1(hiddenSize, visibleSize);
	W1 = MReshape(M_optTheta.MSub(0, 0, hiddenSize * visibleSize, 1), hiddenSize, visibleSize);
	MWriteFile("W1g.dat", W1);

	//cout << "cost = " << cost << endl;
	return 0;
}
