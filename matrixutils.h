#pragma once
#include "def.h"

#define m_abs(a) (((a) > 0.0 ? (a) : (-(a))))
#define m_max(a, b) (((a) > (b) ? (a) : (b)))
#define m_min(a, b) (((a) < (b) ? (a) : (b)))

double RandomMatrixValue ();
double **AllocateMatrix (int rows, int cols);
double **CopyMatrix (double **matrix, int rows, int cols);
double **TransformMatrixByRowOrder (double **matrix, int rows, int cols, int *rowOrder);
double **TransformMatrixByColOrder (double **matrix, int rows, int cols, int *colOrder);
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
void AddMultipliedRowPart (double **matrix, int rows, int cols, int dst, int src, double modifier, int startFrom);
void AddMultipliedCol (double **matrix, int rows, int cols, int dst, int src, double modifier);
