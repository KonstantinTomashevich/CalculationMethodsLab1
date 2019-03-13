#include "lup.h"
#include "matrixutils.h"

#include <stdlib.h>
#include <stdio.h>

static void SelectMax (double **A, int matrixSize, int *PrOrder, int *PcOrder, int step)
{
    double max = A[step][step];
    int maxRow = step;
    int maxCol = step;

    for (int row = step; row < matrixSize; ++row)
    {
        for (int col = step; col < matrixSize; ++col)
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
        SwapRows (A, matrixSize, matrixSize, step, maxRow);
        int temp = PrOrder[step];
        PrOrder[step] = PrOrder[maxRow];
        PrOrder[maxRow] = temp;
    }

    if (maxCol != step)
    {
        SwapCols (A, matrixSize, matrixSize, step, maxCol);
        int temp = PcOrder[step];
        PcOrder[step] = PcOrder[maxCol];
        PcOrder[maxCol] = temp;
    }
}

bool BuildLUP (double **A, int matrixSize, int *PrOrder, int *PcOrder)
{
    for (int index = 0; index < matrixSize; ++index)
    {
        PrOrder [index] = index;
        PcOrder [index] = index;
    }

    for (int step = 0; step < matrixSize; ++step)
    {
        SelectMax (A, matrixSize, PrOrder, PcOrder, step);
        if (A[step][step] == 0.0 || A[step][step] == -0.0)
        {
            return false;
        }

        if (step != matrixSize - 1)
        {
            for (int row = step + 1; row < matrixSize; ++row)
            {
                double modifier = -A[row][step] / A[step][step];
                AddMultipliedRowPart (A, matrixSize, matrixSize, row, step, modifier, step);
                A[row][step] = -modifier;
            }
        }
    }

    return true;
}

bool SolveLUP (double **LU, double **B, int matrixSize, int results, int *PrOrder, int *PcOrder, double ***X)
{
    // Calculate Y from L*Y = Pr*B.
    double **currentResult = TransformMatrixByRowOrder (B, matrixSize, results, PrOrder);

    for (int step = 0; step < matrixSize - 1; ++step)
    {
        for (int row = step + 1; row < matrixSize; ++row)
        {
            AddMultipliedRow (currentResult, matrixSize, results, row, step, -LU[row][step]);
        }
    }

    // Calculate Z from U*Z = Y.
    for (int step = matrixSize - 1; step >= 0; --step)
    {
        if (step > 0)
        {
            for (int row = step - 1; row >= 0; --row)
            {
                AddMultipliedRow (currentResult, matrixSize, results, row, step, -LU[row][step] / LU[step][step]);
            }
        }

        MultiplyRow (currentResult, matrixSize, results, step, 1.0/LU[step][step]);
    }

    // Calculate X from Pc * X = Z.
    *X = AllocateMatrix (matrixSize, results);
    for (int row = 0; row < matrixSize; ++row)
    {
        for (int col = 0; col < results; ++col)
        {
            (*X)[row][col] = currentResult[PcOrder[row]][col];
        }
    }

    FreeMatrix (currentResult, matrixSize, results);
    return true;
}
