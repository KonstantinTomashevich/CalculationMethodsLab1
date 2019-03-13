#include <stdio.h>
#include <stdlib.h>

#include "matrixutils.h"
#include "mtwister.h"
#include "gaussjordan.h"
#include "gauss.h"
#include "lup.h"

#define RUN_COUNT 100
extern MTRand *GlobalRand;
double gaussTotalMaxDiff = 0.0;
double gaussTotalMinDiff = INFINITY;
double gaussTotalAverageDiff = 0.0;


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
        double minDifference = INFINITY;
        double averageDifference = 0.0;

        for (int index = 0; index < MATRIX_SIZE; ++index)
        {
            double currentDiff = m_abs (X[Xi[index]][0] - copyB[index][0]);
            maxDifference = m_max (maxDifference, currentDiff);
            minDifference = m_min (m_abs (minDifference), m_abs (currentDiff));
            averageDifference += currentDiff / MATRIX_SIZE;
        }

        printf ("Gauss max difference: %20.13lf.\n", maxDifference);
        printf ("Gauss min difference: %20.13lf.\n", minDifference);
        printf ("Gauss average difference: %20.13lf.\n", averageDifference);

        gaussTotalMaxDiff = m_max (gaussTotalMaxDiff, maxDifference);
        gaussTotalMinDiff = m_min (gaussTotalMaxDiff, minDifference);
        gaussTotalAverageDiff += averageDifference / RUN_COUNT;
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

    for (int index = 0; index < RUN_COUNT; ++index)
    {
        printf ("\n### Run %d ### \n", index);
        MainCycle ();
    }

    printf ("\n### Report ###\n");
    printf ("Gauss max difference: %20.13lf.\n", gaussTotalMaxDiff);
    printf ("Gauss min difference: %20.13lf.\n", gaussTotalMinDiff);
    printf ("Gauss average difference: %20.13lf.\n", gaussTotalAverageDiff);

    free (GlobalRand);
    return 0;
}
