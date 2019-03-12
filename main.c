#include <stdio.h>
#include <stdlib.h>

#include "matrixutils.h"
#include "mtwister.h"
#include "gaussjordan.h"
#include "gauss.h"

extern MTRand *GlobalRand;
void CalculateConditionNumber (double **A)
{
    double **copyA = CopyMatrix (A, MATRIX_SIZE, MATRIX_SIZE);
    double **antiA = AllocateMatrix (MATRIX_SIZE, MATRIX_SIZE);

    if (!GaussJordanAlgo (copyA, antiA, MATRIX_SIZE, MATRIX_SIZE))
    {
        printf ("Unable to calculate condition number!\n");
    }
    else
    {
        double aMax = m_abs (A[0][0]);
        double antiAMax = m_abs (antiA[0][0]);

        for (int row = 0; row < MATRIX_SIZE; ++row)
        {
            for (int col = 0; col < MATRIX_SIZE; ++col)
            {
                aMax = m_max (aMax, m_abs (A[row][col]));
                antiAMax = m_max (antiAMax, m_abs (antiA[row][col]));
            }
        }

        printf ("Condition number: %lf.\n", aMax * antiAMax);
    }

    FreeMatrix (copyA, MATRIX_SIZE, MATRIX_SIZE);
    FreeMatrix (antiA, MATRIX_SIZE, MATRIX_SIZE);
}

void FindGaussSolutionAndPrintDiff (double **A, double **B, double **X)
{
    double **copyA = CopyMatrix (A, MATRIX_SIZE, MATRIX_SIZE);
    double **copyB = CopyMatrix (B, MATRIX_SIZE, 1);
    int *Xi = calloc (MATRIX_SIZE, sizeof (int));

    if (!Gauss (copyA, copyB, Xi, MATRIX_SIZE, MATRIX_SIZE, 1))
    {
        printf ("Unable to solve system!\n");
    }
    else
    {
        double maxDifference = 0.0;
        double averageDifference = 0.0;

        for (int index = 0; index < MATRIX_SIZE; ++index)
        {
            double currentDiff = m_abs (X[Xi[index]][0] - copyB[index][0]);
            maxDifference = m_max (maxDifference, currentDiff);
            averageDifference += currentDiff / MATRIX_SIZE;
            printf ("X%d diff (X/Gauss): %20.13lf %20.13lf\n", Xi[index], X[Xi[index]][0], copyB[index][0]);
        }

        printf ("Gauss max difference: %20.13lf.\n", maxDifference);
        printf ("Gauss average difference: %20.13lf.\n", averageDifference);
    }

    FreeMatrix (copyA, MATRIX_SIZE, MATRIX_SIZE);
    FreeMatrix (copyB, MATRIX_SIZE, 1);
    free (Xi);
}

void MainCycle ()
{
    double **A = AllocateMatrix (MATRIX_SIZE, MATRIX_SIZE);
    FillTaskSpecificMatrix (A);

    double **X = AllocateMatrix (MATRIX_SIZE, 1);
    FillDefaultMatrix (X, MATRIX_SIZE, 1);

    double **B = AllocateMatrix (MATRIX_SIZE, 1);
    MultiplyMatrices (A, X, B, MATRIX_SIZE, MATRIX_SIZE, 1);

    CalculateConditionNumber (A);
    FindGaussSolutionAndPrintDiff (A, B, X);

    FreeMatrix (A, MATRIX_SIZE, MATRIX_SIZE);
    FreeMatrix (X, MATRIX_SIZE, 1);
    FreeMatrix (B, MATRIX_SIZE, 1);
}

int main ()
{
    GlobalRand = malloc (sizeof (MTRand));
    *GlobalRand = SeedRand (1377);

    MainCycle ();
    free (GlobalRand);
    return 0;
}
