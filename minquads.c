#include "minquads.h"
#include "matrixutils.h"
#include "householder.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

bool SolveMinQuads (double **A, int rows, int cols, double **B, int results, double ***X)
{
    double **At = TransposeMatrix (A, rows, cols);
    double **AtxA = AllocateMatrix (cols, cols);
    *X = AllocateMatrix (cols, results);

    MultiplyMatrices (At, A, AtxA, cols, rows, cols);
    MultiplyMatrices (At, B, *X, cols, rows, results);
    SolveHouseholder (AtxA, cols, *X, results);

    FreeMatrix (At, cols, rows);
    FreeMatrix (AtxA, cols, cols);
    return true;
}
