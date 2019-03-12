#include "gaussjordan.h"
#include "matrixutils.h"

static bool FixA (double **A, double **antiA, int rows, int cols, int step)
{
    if (A [step] [step] == 0)
    {
        if (step == rows - 1)
        {
            return false;
        }

        int toChange = -1;
        for (int row = step + 1; row < rows; ++row)
        {
            if (A [row] [step] != 0)
            {
                toChange = row;
                break;
            }
        }

        if (toChange == -1)
        {
            return false;
        }

        SwapRows (A, rows, cols, step, toChange);
        SwapRows (antiA, rows, cols, step, toChange);
    }

    return true;
}

bool GaussJordanAlgo (double **A, double **antiA, int rows, int cols)
{
    if (rows != cols)
    {
        return false;
    }

    for (int row = 0; row < rows; ++row)
    {
        for (int col = 0; col < cols; ++col)
        {
            antiA [row] [col] = row == col ? 1 : 0;
        }
    }

    for (int step = 0; step < rows; ++step)
    {
        if (!FixA (A, antiA, rows, cols, step))
        {
            return false;
        }

        if (step != rows - 1)
        {
            for (int row = step + 1; row < rows; ++row)
            {
                double modifier = -A [row] [step] / A [step] [step];
                AddMultipliedRow (A, rows, cols, row, step, modifier);
                AddMultipliedRow (antiA, rows, cols, row, step, modifier);
            }
        }
    }

    for (int step = 0; step < rows; ++step)
    {
        double modifier = 1.0 / A [step] [step];
        MultiplyRow (A, rows, cols, step, modifier);
        MultiplyRow (antiA, rows, cols, step, modifier);
    }

    for (int step = rows - 1; step > 0; --step)
    {
        for (int row = step - 1; row >= 0; --row)
        {
            double modifier = -A [row] [step] / A [step] [step];
            AddMultipliedRow (A, rows, cols, row, step, modifier);
            AddMultipliedRow (antiA, rows, cols, row, step, modifier);
        }
    }

    return true;
}
