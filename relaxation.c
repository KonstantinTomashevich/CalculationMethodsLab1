#include "relaxation.h"
#include "matrixutils.h"

#include <stdio.h>
#include <math.h>

#define RELAXATION_SOLVED_COEFF 0.00000000000001

bool SolveRelaxation (double **A, int matrixSize, double **B, int results, double ***X)
{
    // x(n+1) = C * x(n) + B.
    // Transform A to C.
    for (int step = 0; step < matrixSize; ++step)
    {
        if (A[step][step] == 0)
        {
            return false;
        }

        double multiplier = 1.0 / A[step][step];
        MultiplyRow (A, matrixSize, matrixSize, step, -multiplier);
        MultiplyRow (B, matrixSize, results, step, multiplier);
        A[step][step] += 1.0;
    }

    *X = CopyMatrix (B, matrixSize, results);
    double **previousX = AllocateMatrix (matrixSize, matrixSize);
    bool moreStepsNeeded = true;
    const double coeff = (N + 1) / 6.0;

    while (moreStepsNeeded)
    {
        for (int row = 0; row < matrixSize; ++row)
        {
            for (int col = 0; col < results; ++col)
            {
                previousX[row][col] = (*X)[row][col];
            }
        }

        for (int step = 0; step < matrixSize; ++step)
        {
            for (int col = 0; col < results; ++col)
            {
                double value = 0.0;
                for (int index = 0; index < matrixSize; ++index)
                {
                    value += A[step][index] * (*X)[index][col];
                }

                (*X)[step][col] = value + B[step][col];
            }
        }

        double maxDiff = 0.0;
        for (int row = 0; row < matrixSize; ++row)
        {
            for (int col = 0; col < results; ++col)
            {
                (*X)[row][col] = (1 - coeff) * previousX[row][col] + coeff * (*X)[row][col];
                double diff = fabs ((*X)[row][col] - previousX[row][col]);

                if (diff > maxDiff)
                {
                    maxDiff = diff;
                }
            }
        }

        moreStepsNeeded = maxDiff > RELAXATION_SOLVED_COEFF;
    }

    return true;
}
