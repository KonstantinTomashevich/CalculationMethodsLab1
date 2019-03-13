#pragma once
#include <stdbool.h>

/// Result: Diagonal and upper triangle of A is U, lower triangle is L (L diagonal is ones).
bool BuildLUP (double **A, int matrixSize, int *PrOrder, int *PcOrder);
/// LU, PrOrder and PcOrder must NOT be changed during execution!
/// B wouldn't be changed.
/// X is output pointer to created solution matrix.
bool SolveLUP (double **LU, double **B, int matrixSize, int results, int *PrOrder, int *PcOrder, double ***X);
