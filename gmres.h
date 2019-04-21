#pragma once
#include <stdbool.h>

/// Note: B should be [matrixSize]x1 sized matrix.
bool SolveGMRES (double **A, int matrixSize, double **B, double ***X);
