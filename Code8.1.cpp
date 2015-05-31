#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "mat.h"

using namespace std;

const int patchSize = 12;
const int numPatches = 10000;

Mat sampleIMAGES()
{
	int imageCount = 10;

	int row = 512;
	int col = 512;

	double** IMAGES = MReadFile("IMAGERAW.dat", imageCount * row, col);
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
	M_patches = M_patches - MRepmat(MMean(M_patches), M_patches.row, 1);

	return M_patches;
}

int main()
{
	//double** patches = MReadFile("patches.dat", patchSize * patchSize, numPatches);
	Mat M_patches(patchSize * patchSize, numPatches);
	M_patches = sampleIMAGES();
	//MWriteFile("patches.dat", M_patches);

	Mat M_sigma(patchSize * patchSize, patchSize * patchSize);
	M_sigma = (1.0 / numPatches) * M_patches * MTrans(M_patches);

	Mat M_U(patchSize * patchSize, patchSize * patchSize);
	Mat M_S(patchSize * patchSize, patchSize * patchSize);
	Mat M_V(patchSize * patchSize, patchSize * patchSize);

	Mat M_A(patchSize * patchSize, patchSize * patchSize);

	M_A = M_sigma.SVD(M_U, M_S, M_V);

	Mat M_xRot(patchSize * patchSize, numPatches);
	M_xRot = MTrans(M_U) * M_patches;

	Mat M_covar(patchSize * patchSize, patchSize * patchSize);
	M_covar = (1.0 / numPatches) * M_xRot * MTrans(M_xRot);

	Mat M_ss(patchSize * patchSize, 1);
	M_ss = M_covar.Diag();

	Mat k(1, 1);
	k = MSum(MCusum(M_ss) / MSum(M_ss) <= 0.99);

	Mat M_Un(patchSize * patchSize, numPatches);
	M_Un.MPub(MTrans(M_U.MSub(0, 0, patchSize * patchSize, k.data[0][0])) * M_patches, 0, 0);
	
	Mat M_xHat(patchSize * patchSize, numPatches);
	M_xHat = M_U * M_Un;

	double epsilon = 0.1;
	Mat M_xPCAWhite(patchSize * patchSize, numPatches);
	M_xPCAWhite = MDiag(1.0 / MDPow(M_A.Diag() + epsilon, 0.5)) * MTrans(M_U) * M_patches;

	M_covar = (1.0 / numPatches) * M_xPCAWhite * MTrans(M_xPCAWhite);

	Mat M_xZCAWhite(patchSize * patchSize, numPatches);
	M_xZCAWhite = M_U * M_xPCAWhite;

	MWriteFile("xZCAWhite.dat", M_xZCAWhite);

	return 0;
}
