#include "gmres.h"
#include "matrixutils.h"
#include "minquads.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define EPSILON 0.00001
bool SolveGMRES (double **A, int matrixSize, double **B, double ***X)
{
    double **K = AllocateMatrix (matrixSize, matrixSize);
    double **AK = AllocateMatrix (matrixSize, matrixSize);
    double **C;
    double **previousX = AllocateMatrix (matrixSize, 1);
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

        // TODO: Check diff AX and B instead.
        if (iteration > 0)
        {
            double diff = 0.0;
            for (int row = 0; row < matrixSize; ++row)
            {
                diff = fmax (diff, fabs (previousX[row][0] - (*X)[row][0]));
            }

            if (diff < EPSILON)
            {
                solved = true;
            }
        }

        CopyMatrixInto (*X, matrixSize, 1, previousX);
        FreeMatrix (C, iteration + 1, 1);
        ++iteration;
    }

    FreeMatrix (K, matrixSize, matrixSize);
    FreeMatrix (AK, matrixSize, matrixSize);
    FreeMatrix (previousX, matrixSize, 1);
    return solved;
}
