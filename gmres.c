#include "gmres.h"
#include "matrixutils.h"
#include "minquads.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define EPSILON 0.0001
bool SolveGMRES (double **A, int matrixSize, double **B, double ***X)
{
    double **K = AllocateMatrix (matrixSize, matrixSize);
    double **AK = AllocateMatrix (matrixSize, matrixSize);
    double **C;
    double **sample = AllocateMatrix (matrixSize, 1);
    *X = AllocateMatrix (matrixSize, 1);

    bool solved = false;
    int iteration = 0;

    while (iteration < matrixSize - 2 && !solved)
    {
        if (iteration == 0)
        {
            for (int row = 0; row < matrixSize; ++row)
            {
                K[row][0] = B[row][0];
            }
        }
        else
        {
            for (int row = 0; row < matrixSize; ++row)
            {
                K[row][iteration] = AK[row][iteration - 1];
            }
        }

        for (int row = 0; row < matrixSize; ++row)
        {
            double value = 0.0;
            for (int mulIndex = 0; mulIndex < matrixSize; ++mulIndex)
            {
                value += A[row][mulIndex] * K[mulIndex][iteration];
            }

            AK[row][iteration] = value;
        }

        SolveMinQuads (AK, matrixSize, iteration + 1, B, 1, &C);
        MultiplyMatrices (K, C, *X, matrixSize, iteration + 1, 1);
        MultiplyMatrices (A, *X, sample, matrixSize, matrixSize, 1);


        double normal = 0.0;
        for (int row = 0; row < matrixSize; ++row)
        {
            normal = pow (sample[row][0] - B[row][0], 2.0);
        }

        normal = sqrt (normal);
        if (normal < EPSILON)
        {
            solved = true;
        }

        FreeMatrix (C, iteration + 1, 1);
        ++iteration;
    }

    FreeMatrix (K, matrixSize, matrixSize);
    FreeMatrix (AK, matrixSize, matrixSize);
    FreeMatrix (sample, matrixSize, 1);
    return solved;
}
