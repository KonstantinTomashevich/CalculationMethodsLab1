#include "minquads.h"
#include "matrixutils.h"
#include "gaussjordan.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

bool SolveMinQuads (double **A, int rows, int cols, double **B, int results, double ***X)
{
    double **At = TransposeMatrix (A, rows, cols);
    double **AtxA = AllocateMatrix (cols, cols);
    double **antiAtxA = AllocateMatrix (cols, cols);
    double **thirdStep = AllocateMatrix (cols, rows); // antiAtxA * At
    *X = AllocateMatrix (cols, results);

    MultiplyMatrices (At, A, AtxA, cols, rows, cols);
    GaussJordanAlgo (AtxA, antiAtxA, cols, cols);
    MultiplyMatrices (antiAtxA, At, thirdStep, cols, cols, rows);
    MultiplyMatrices (thirdStep, B, *X, cols, rows, results);

    FreeMatrix (At, cols, rows);
    FreeMatrix (AtxA, cols, cols);
    FreeMatrix (antiAtxA, cols, cols);
    FreeMatrix (thirdStep, cols, cols);
    return true;
}
