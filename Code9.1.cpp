#include <iostream>
#include <cstdlib>
#include <ctime>
#include "mat.h"
#include "minFunc.h"

using namespace std;

const int inputSize = 28 * 28;
const int numClasses = 10;

const int numCases = 60000;

const double lambda = 1e-4;
const double sparsityParam = 0.01;
const double beta = 3;

double softmaxCost
(Mat& theta, 
 const int& numClasses, 
 const int& inputSize, 
 const double& lambda,
 const Mat& data,
 const Mat& labels, 
 Mat& grad)
{
	Mat x(numClasses, inputSize);
	x = MReshape(theta, numClasses, inputSize);

	Mat groundTruth(numClasses, numCases);
	for(int i = 0; i < numCases; i++)
	{
		double idx = labels.data[i][0];
		int index = int(idx);
		groundTruth[index][i] = 1;
	}

	Mat cost(1, 1);
	Mat M(numClasses, numCases);
	M = x * data;

	M -= MRepmat(Max(M), numClasses, 1);
	
	M = MExp(M);
	M /= MRepmat(MSum(M), numClasses, 1);

	//cost = -1.0 / numCases * MTrans(groundTruth.All()) * MLog(M.All())
	//	+ lambda / 2.0 * MSum(MDPow(x.All(), 2));

	Mat thetagrad(numClasses, inputSize);
	thetagrad = -1.0 / numCases * (groundTruth - M) * MTrans(data) + lambda * x;
	grad = thetagrad.All();

	return cost[0][0];
}

Mat softmaxTrain
(const int& inputSize, 
 const int& numClasses, 
 const double& lambda, 
 const Mat& data,
 const Mat& labels
 )
{
	double** th = MReadFile("small_theta_init01.dat", numClasses * inputSize, 1);
	Mat theta(numClasses * inputSize, 1, th);
	/*for(int i = 0; i < numClasses * inputSize; i++)
	{
		theta[i][0] = double(rand() % 1000) / 1000 * 0.005;
	}*/

	Mat M_grad(numClasses * inputSize, 1);

	softmaxCost(theta, numClasses, inputSize, lambda, data, labels, M_grad);

	Mat M_opttheta(numClasses * inputSize, 1);

	minFuncOp opt;
	opt.maxIter = 100;
	opt.type = 2;

	M_opttheta = minFunc(&softmaxCost, theta, numClasses, inputSize, lambda, data, labels, M_grad, opt);
	return MReshape(M_opttheta, numClasses, inputSize);
}

Mat softmaxPredict(const Mat& M_model, const Mat& data)
{
	Mat M_order(1, data.col);
	Max(M_model * data, M_order);
	return M_order;
}

int main()
{
	double** trainImages = MReadFile("trainImages.dat", inputSize, numCases);
	double** trainLabels = MReadFile("trainLabels.dat", numCases, 1);
	Mat M_train(inputSize, numCases, trainImages);

	Mat M_label(numCases, 1, trainLabels);

	Mat M_model(numClasses, inputSize);
	M_model = softmaxTrain(inputSize, numClasses, lambda, M_train, M_label);

	Mat M_res(1, 10000);

	double** testImages = MReadFile("testImages.dat", inputSize, 10000);
	Mat M_test(inputSize, 10000, testImages);
	double** testLabels = MReadFile("testLabels.dat", 10000, 1);
	Mat M_testl(10000, 1, testLabels);

	M_res = softmaxPredict(M_model, M_test);
	Mat M_acc(1, 1);
	M_acc = (MSum(M_res.All() == M_testl)).All() / 10000;
	M_acc.MShow();

	return 0;
}
