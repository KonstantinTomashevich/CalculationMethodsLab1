#include "gauss.h"
#include "matrixutils.h"

#include <stdlib.h>
#include <stdio.h>

static void SelectMax (double **A, double **B, int *Xi, int rows, int cols, int results, int step)
{
    double max = A[step][step];
    int maxRow = step;
    int maxCol = step;

    for (int row = step; row < rows; ++row)
    {
        for (int col = step; col < cols; ++col)
        {
            if (m_abs (A[row][col]) > m_abs (max))
            {
                max = A[row][col];
                maxRow = row;
                maxCol = col;
            }
        }
    }

    if (maxRow != step)
    {
        SwapRows (A, rows, cols, step, maxRow);
        SwapRows (B, rows, results, step, maxRow);
    }

    if (maxCol != step)
    {
        SwapCols (A, rows, cols, step, maxCol);
        int temp = Xi[step];
        Xi[step] = Xi[maxCol];
        Xi[maxCol] = temp;
    }
}

bool Gauss (double **A, double **B, int *Xi, int rows, int cols, int results)
{
    for (int index = 0; index < cols; ++index)
    {
        Xi[index] = index;
    }

    for (int step = 0; step < rows; ++step)
    {
        SelectMax (A, B, Xi, rows, cols, results, step);
        if (step != rows - 1 && A[step][step] != 0.0 && A[step][step] != -0.0)
        {
            for (int row = step + 1; row < rows; ++row)
            {
                double modifier = -A[row][step] / A[step][step];
                AddMultipliedRow (A, rows, cols, row, step, modifier);
                AddMultipliedRow (B, rows, results, row, step, modifier);
            }
        }
    }

    for (int step = rows - 1; step > 0; --step)
    {
        if (A[step][step] != 0.0 && A[step][step] != -0.0)
        {
            for (int row = step - 1; row >= 0; --row)
            {
                double modifier = -A[row][step] / A[step][step];
                A[row][step] += modifier * A[step][step];
                AddMultipliedRow (B, rows, results, row, step, modifier);
            }
        }
    }

    int firstZero = -1;
    for (int step = 0; step < rows; ++step)
    {
        if (A[step][step] != 0.0 && A[step][step] != -0.0)
        {
            double modifier = 1.0 / A[step][step];
            MultiplyRow (A, rows, cols, step, modifier);
            MultiplyRow (B, rows, results, step, modifier);
        }
        else if (firstZero == -1)
        {
            firstZero = step;
        }
    }

    if (firstZero != -1)
    {
        for (int row = firstZero; row < rows; ++row)
        {
            for (int col = 0; col < results; ++col)
            {
                if (B[row][col] != 0.0 || B[row][col] != -0.0)
                {
                    return false;
                }
            }
        }
    }

    return true;
}

void PrintGaussSolution (double **sA, double **sB, int *Xi, int rows, int cols, int results)
{
    int *positions = calloc (cols, sizeof (int));
    for (int index = 0; index < cols; ++index)
    {
        positions[Xi[index]] = index;
    }

    for (int index = 0; index < cols; ++index)
    {
        int position = positions[index];
        printf ("X%4d: ", index);

        if (sA[position][position] == 0.0 || sA[position][position] == -0.0 || index > rows)
        {
            printf ("%s", "ERROR: found free variable, not supported yet!");
        }
        else
        {
            for (int col = 0; col < results; ++col)
            {
                printf ("%20lf ", sB[position][col]);
            }
        }

        printf ("\n");
    }
}
