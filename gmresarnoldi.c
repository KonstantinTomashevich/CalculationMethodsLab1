#include "gmresarnoldi.h"
#include "matrixutils.h"
#include "minquads.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define EPSILON 0.00001
bool SolveGMRESArnoldi (double **A, int matrixSize, double **B, double ***X)
{
    double **Q = AllocateMatrix (matrixSize, matrixSize);
    double **H = AllocateMatrix (matrixSize, matrixSize);
    double **minquadsB = AllocateMatrix (matrixSize, 1);
    double **C;
    double **previousX = AllocateMatrix (matrixSize, 1);
    *X = AllocateMatrix (matrixSize, 1);
    bool solved = false;

    for (int iteration = 0; iteration < matrixSize; ++iteration)
    {
        if (iteration == 0)
        {
            double normal = ColumnQuadricNormal (B, matrixSize, 0);
            for (int row = 0; row < matrixSize; ++row)
            {
                Q[row][0] = B[row][0] / normal;
            }

            minquadsB[0][0] = normal;
            continue;
        }

        for (int row = 0; row < matrixSize; ++row)
        {
            double value = 0.0;
            for (int mulIndex = 0; mulIndex < matrixSize; ++mulIndex)
            {
                value += A[row][mulIndex] * Q[mulIndex][iteration - 1];
            }

            Q[row][iteration] = value;
        }

        for (int normalizationStep = 0; normalizationStep < iteration; ++normalizationStep)
        {
            H[normalizationStep][iteration - 1] =
                ColumnsScalarMultiplication (Q, matrixSize, iteration, normalizationStep);

            for (int row = 0; row < matrixSize; ++row)
            {
                Q[row][iteration] = Q[row][iteration] -
                    H[normalizationStep][iteration - 1] * Q[row][normalizationStep];
            }
        }

        H[iteration][iteration - 1] = ColumnQuadricNormal (Q, matrixSize, iteration);
        if (H[iteration][iteration - 1] == 0)
        {
            break;
        }

        for (int row = 0; row < matrixSize; ++row)
        {
            Q[row][iteration] /= H[iteration][iteration - 1];
        }

        SolveMinQuads (H, iteration + 1, iteration, minquadsB, 1, &C);
        MultiplyMatrices (Q, C, *X, matrixSize, iteration, 1);

        if (iteration > 1)
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
        FreeMatrix (C, iteration, 1);

        if (solved)
        {
            break;
        }
    }

    FreeMatrix (Q, matrixSize, matrixSize);
    FreeMatrix (H, matrixSize, matrixSize);
    FreeMatrix (previousX, matrixSize, 1);
    return solved;
}
