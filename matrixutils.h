#pragma once
#include "def.h"

double RandomMatrixValue ();
double **AllocateMatrix (int rows, int cols);
double **CopyMatrix (double **matrix, int rows, int cols);
void FreeMatrix (double **matrix, int rows, int cols);

void FillTaskSpecificMatrix (double **matrix);
void FillDefaultMatrix (double **matrix, int rows, int cols);
void PrintMatrix (double **matrix, int rows, int cols);

void MultiplyRow (double **matrix, int rows, int cols, int row, double by);
void MultiplyCol (double **matrix, int rows, int cols, int col, double by);
void SwapRows (double **matrix, int rows, int cols, int row1, int row2);
void SwapCols (double **matrix, int rows, int cols, int col1, int col2);

void MultiplyMatrices (double **first, double **second, double **output, int firstRows, int firstCols, int secondsCols);
void AddMultipliedRow (double **matrix, int rows, int cols, int dst, int src, double modifier);
void AddMultipliedCol (double **matrix, int rows, int cols, int dst, int src, double modifier);