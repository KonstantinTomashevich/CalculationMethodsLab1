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
