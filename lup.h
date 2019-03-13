#pragma once
#include <stdbool.h>

/// Result: Diagonal and upper triangle of A is U, lower triangle is L (L diagonal is ones).
bool BuildLUP (double **A, int matrixSize, int *PrOrder, int *PcOrder);
