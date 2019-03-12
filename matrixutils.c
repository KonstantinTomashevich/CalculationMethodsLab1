#include "matrixutils.h"
#include "mtwister.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

extern MTRand *GlobalRand;
static double MinValue = 0.0;
static double MaxValue = 0.0;
static double Addition = 0.0;
static double Numberer = 0.0;

double RandomMatrixValue ()
{
    if (MinValue == MaxValue)
    {
        MaxValue = pow (2.0, N / 4.0);
        MinValue = -MaxValue;
        Addition = pow (10.0, -13.0);
        Numberer = pow (10.0, 13.0);
    }

    double result = GenRand (GlobalRand) * (MaxValue - MinValue) + MinValue;
    long long asNumber = (long long) (result * Numberer);

    if (asNumber % 10 == 0)
    {
        result += Addition;
    }

    return result;
}

double **AllocateMatrix (int rows, int cols)
{
    double **matrix = calloc (rows, sizeof (double *));
    for (int index = 0; index < rows; ++index)
    {
        matrix[index] = calloc (cols, sizeof (double));
    }

    return matrix;
}

double **CopyMatrix (double **matrix, int rows, int cols)
{
    double **copy = calloc (rows, sizeof (double *));
    for (int row = 0; row < rows; ++row)
    {
        copy[row] = calloc (cols, sizeof (double));
        for (int col = 0; col < cols; ++col)
        {
            copy[row][col] = matrix[row][col];
        }
    }

    return copy;
}

void FreeMatrix (double **matrix, int rows, int cols)
{
    for (int row = 0; row < rows; ++row)
    {
        free (matrix[row]);
    }

    free (matrix);
}

void FillTaskSpecificMatrix (double **matrix)
{
    for (int col = 1; col < MATRIX_SIZE; ++col)
    {
        for (int row = 0; row < col; ++row)
        {
            double value = RandomMatrixValue ();
            matrix[row][col] = value;
            matrix[col][row] = value;
        }
    }

    for (int diag = 0; diag < MATRIX_SIZE; ++diag)
    {
        double value = 0.0;
        for (int col = 0; col < MATRIX_SIZE; ++col)
        {
            if (col != diag)
            {
                value += matrix[diag][col];
            }
        }

        matrix[diag][diag] = value;
    }
}

void FillDefaultMatrix (double **matrix, int rows, int cols)
{
    for (int row = 0; row < rows; ++row)
    {
        for (int col = 0; col < cols; ++col)
        {
            matrix[row][col] = RandomMatrixValue ();
        }
    }
}

void PrintMatrix (double **matrix, int rows, int cols)
{
    for (int row = 0; row < rows; ++row)
    {
        for (int col = 0; col < cols; ++col)
        {
            printf ("%20.13lf ", matrix[row][col]);
        }

        printf ("\n");
    }
}

void MultiplyRow (double **matrix, int rows, int cols, int row, double by)
{
    for (int col = 0; col < cols; ++col)
    {
        matrix[row][col] *= by;
    }
}

void MultiplyCol (double **matrix, int rows, int cols, int col, double by)
{
    for (int row = 0; row < rows; ++row)
    {
        matrix[row][col] *= by;
    }
}

void SwapRows (double **matrix, int rows, int cols, int row1, int row2)
{
    double *temp = matrix[row1];
    matrix[row1] = matrix[row2];
    matrix[row2] = temp;
}

void SwapCols (double **matrix, int rows, int cols, int col1, int col2)
{
    for (int row = 0; row < rows; ++row)
    {
        double temp = matrix[row][col1];
        matrix[row][col1] = matrix[row][col2];
        matrix[row][col2] = temp;
    }
}

void MultiplyMatrices (double **first, double **second, double **output, int firstRows, int firstCols, int secondsCols)
{
    for (int row = 0; row < firstRows; ++row)
    {
        for (int col = 0; col < secondsCols; ++col)
        {
            double value = 0.0;
            for (int iterator = 0; iterator < firstCols; ++iterator)
            {
                value += first[row][iterator] * second[iterator][col];
            }

            output[row][col] = value;
        }
    }
}

void AddMultipliedRow (double **matrix, int rows, int cols, int dst, int src, double modifier)
{
    for (int col = 0; col < cols; ++col)
    {
        matrix[dst][col] += matrix[src][col] * modifier;
    }
}

void AddMultipliedCol (double **matrix, int rows, int cols, int dst, int src, double modifier)
{
    for (int row = 0; row < rows; ++row)
    {
        matrix[row][dst] += matrix[row][src] * modifier;
    }
}
