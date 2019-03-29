#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "matrixutils.h"
#include "mtwister.h"
#include "gaussjordan.h"
#include "gauss.h"
#include "lup.h"
#include "cholesky.h"
#include "relaxation.h"

#define RUN_COUNT 100
extern MTRand *GlobalRand;
double minConditionNumber = INFINITY;
double maxConditionNumber = 0.0;
double averageMatrixElement = 0.0;

clock_t totalGaussJordan = 0;

double gaussTotalMaxDiff = 0.0;
double gaussTotalMinDiff = INFINITY;
double gaussTotalAverageDiff = 0.0;

clock_t totalGauss = 0;

clock_t totalLUPBuild = 0;

double lupTotalMaxDiff = 0.0;
double lupTotalMinDiff = INFINITY;
double lupTotalAverageDiff = 0.0;

clock_t totalLUPSolve = 0;

double choleskyTotalMaxDiff = 0.0;
double choleskyTotalMinDiff = INFINITY;
double choleskyTotalAverageDiff = 0.0;

clock_t totalCholeskySolve = 0;

double relaxationTotalMaxDiff = 0.0;
double relaxationTotalMinDiff = INFINITY;
double relaxationTotalAverageDiff = 0.0;

clock_t totalRelaxation = 0;

void CalculateConditionNumber (double **A)
{
    double **copyA = CopyMatrix (A, MATRIX_SIZE, MATRIX_SIZE);
    double **antiA = AllocateMatrix (MATRIX_SIZE, MATRIX_SIZE);
    clock_t begin = clock ();

    if (!GaussJordanAlgo (copyA, antiA, MATRIX_SIZE, MATRIX_SIZE))
    {
        printf ("Unable to calculate condition number!\n");
    }
    else
    {
        totalGaussJordan += clock () - begin;
        double aMax = m_abs (A[0][0]);
        double antiAMax = m_abs (antiA[0][0]);
        double average = 0.0;

        for (int row = 0; row < MATRIX_SIZE; ++row)
        {
            for (int col = 0; col < MATRIX_SIZE; ++col)
            {
                average += A[row][col] / (MATRIX_SIZE * MATRIX_SIZE);
                aMax = m_max (aMax, m_abs (A[row][col]));
                antiAMax = m_max (antiAMax, m_abs (antiA[row][col]));
            }
        }

        minConditionNumber = m_min (minConditionNumber, aMax * antiAMax);
        maxConditionNumber = m_max (maxConditionNumber, aMax * antiAMax);
        averageMatrixElement += average / RUN_COUNT;
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
    clock_t begin = clock ();

    if (!Gauss (copyA, copyB, Xi, MATRIX_SIZE, MATRIX_SIZE, 1))
    {
        printf ("Unable to solve system!\n");
    }
    else
    {
        totalGauss += clock () - begin;
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

        printf ("Gauss max difference: %23.16lf.\n", maxDifference);
        printf ("Gauss min difference: %23.16lf.\n", minDifference);
        printf ("Gauss average difference: %23.16lf.\n", averageDifference);

        gaussTotalMaxDiff = m_max (gaussTotalMaxDiff, maxDifference);
        gaussTotalMinDiff = m_min (gaussTotalMaxDiff, minDifference);
        gaussTotalAverageDiff += averageDifference / RUN_COUNT;
    }

    FreeMatrix (copyA, MATRIX_SIZE, MATRIX_SIZE);
    FreeMatrix (copyB, MATRIX_SIZE, 1);
    free (Xi);
}

void FindLUPSolutionAndPrintDiff (double **A, double **B, double **X)
{
    double **LU = CopyMatrix (A, MATRIX_SIZE, MATRIX_SIZE);
    double **copyB = CopyMatrix (B, MATRIX_SIZE, 1);
    int PrOrder [MATRIX_SIZE];
    int PcOrder [MATRIX_SIZE];
    clock_t begin = clock ();

    if (!BuildLUP (LU, MATRIX_SIZE, PrOrder, PcOrder))
    {
        printf ("Unable to build LUP!\n");
    }
    else
    {
        totalLUPBuild += clock () - begin;
        double **builtX;
        begin = clock ();

        if (!SolveLUP (LU, copyB, MATRIX_SIZE, 1, PrOrder, PcOrder, &builtX))
        {
            printf ("Unable to solve system with LUP!\n");
        }
        else
        {
            totalLUPSolve += clock () - begin;
            double maxDifference = 0.0;
            double minDifference = INFINITY;
            double averageDifference = 0.0;

            for (int index = 0; index < MATRIX_SIZE; ++index)
            {
                double currentDiff = m_abs (X[index][0] - builtX[index][0]);
                maxDifference = m_max (maxDifference, currentDiff);
                minDifference = m_min (m_abs (minDifference), m_abs (currentDiff));
                averageDifference += currentDiff / MATRIX_SIZE;
            }

            printf ("LUP max difference: %23.16lf.\n", maxDifference);
            printf ("LUP min difference: %23.16lf.\n", minDifference);
            printf ("LUP average difference: %23.16lf.\n", averageDifference);

            lupTotalMaxDiff = m_max (lupTotalMaxDiff, maxDifference);
            lupTotalMinDiff = m_min (lupTotalMinDiff, minDifference);
            lupTotalAverageDiff += averageDifference / RUN_COUNT;
            FreeMatrix (builtX, MATRIX_SIZE, MATRIX_SIZE);
        }
    }

    FreeMatrix (LU, MATRIX_SIZE, MATRIX_SIZE);
    FreeMatrix (copyB, MATRIX_SIZE, 1);
}

void FindCholeskySolutionAndPrintDiff (double **A, double **B, double **X)
{
    double **LT = CopyMatrix (A, MATRIX_SIZE, MATRIX_SIZE);
    double **copyB = CopyMatrix (B, MATRIX_SIZE, 1);
    int D [MATRIX_SIZE];
    clock_t begin = clock ();

    if (!BuildCholeskyLT (LT, MATRIX_SIZE, D))
    {
        printf ("Unable to build Cholesky!\n");
    }
    else
    {
        double **builtX;
        if (!SolveCholesky (LT, MATRIX_SIZE, D, copyB, 1, &builtX))
        {
            printf ("Unable to solve system with Cholesky!\n");
        }
        else
        {
            totalCholeskySolve += clock () - begin;
            double maxDifference = 0.0;
            double minDifference = INFINITY;
            double averageDifference = 0.0;

            for (int index = 0; index < MATRIX_SIZE; ++index)
            {
                double currentDiff = m_abs (X[index][0] - builtX[index][0]);
                maxDifference = m_max (maxDifference, currentDiff);
                minDifference = m_min (m_abs (minDifference), m_abs (currentDiff));
                averageDifference += currentDiff / MATRIX_SIZE;
            }

            printf ("Cholesky max difference: %23.16lf.\n", maxDifference);
            printf ("Cholesky min difference: %23.16lf.\n", minDifference);
            printf ("Cholesky average difference: %23.16lf.\n", averageDifference);

            choleskyTotalMaxDiff = m_max (lupTotalMaxDiff, maxDifference);
            choleskyTotalMinDiff = m_min (lupTotalMinDiff, minDifference);
            choleskyTotalAverageDiff += averageDifference / RUN_COUNT;
            FreeMatrix (builtX, MATRIX_SIZE, MATRIX_SIZE);
        }
    }

    FreeMatrix (LT, MATRIX_SIZE, MATRIX_SIZE);
    FreeMatrix (copyB, MATRIX_SIZE, 1);
}

void FindRelaxationSolutionAndPrintDiff (double **A, double **B, double **X)
{
    double **copyA = CopyMatrix (A, MATRIX_SIZE, MATRIX_SIZE);
    double **copyB = CopyMatrix (B, MATRIX_SIZE, 1);
    double **builtX;
    clock_t begin = clock ();

    if (!SolveRelaxation (copyA, MATRIX_SIZE, copyB, 1, &builtX))
    {
        printf ("Unable to solve system!\n");
    }
    else
    {
        totalRelaxation += clock () - begin;
        double maxDifference = 0.0;
        double minDifference = INFINITY;
        double averageDifference = 0.0;

        for (int index = 0; index < MATRIX_SIZE; ++index)
        {
            double currentDiff = m_abs (builtX[index][0] - X[index][0]);
            maxDifference = m_max (maxDifference, currentDiff);
            minDifference = m_min (m_abs (minDifference), m_abs (currentDiff));
            averageDifference += currentDiff / MATRIX_SIZE;
        }

        printf ("Relaxation max difference: %23.16lf.\n", maxDifference);
        printf ("Relaxation min difference: %23.16lf.\n", minDifference);
        printf ("Relaxation average difference: %23.16lf.\n", averageDifference);

        relaxationTotalMaxDiff = m_max (relaxationTotalMaxDiff, maxDifference);
        relaxationTotalMinDiff = m_min (relaxationTotalMaxDiff, minDifference);
        relaxationTotalAverageDiff += averageDifference / RUN_COUNT;
    }

    FreeMatrix (copyA, MATRIX_SIZE, MATRIX_SIZE);
    FreeMatrix (copyB, MATRIX_SIZE, 1);
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
    FindLUPSolutionAndPrintDiff (A, B, X);
    FindCholeskySolutionAndPrintDiff (A, B, X);
    FindRelaxationSolutionAndPrintDiff (A, B, X);

    FreeMatrix (A, MATRIX_SIZE, MATRIX_SIZE);
    FreeMatrix (X, MATRIX_SIZE, 1);
    FreeMatrix (B, MATRIX_SIZE, 1);
}

int main ()
{
    // TODO: After finishing all algos, speed up them by addition of "start from row/col" parameter to util functions.
    GlobalRand = malloc (sizeof (MTRand));
    *GlobalRand = SeedRand (1377);

    for (int index = 0; index < RUN_COUNT; ++index)
    {
        printf ("\n### Run %d ### \n", index);
        MainCycle ();
    }

    printf ("\n### Report ###\n");
    printf ("## 1\nMax condition number: %23.16lf.\n", maxConditionNumber);
    printf ("Min condition number: %23.16lf.\n", minConditionNumber);
    printf ("Average matrix element: %23.16lf.\n\n", averageMatrixElement);

    printf ("## 2\nAverage A^-1 calculation time: %dms.\n\n",
            (int) round (totalGaussJordan * 1000.0 / CLOCKS_PER_SEC / RUN_COUNT));

    printf ("## 3\nGauss max difference: %23.16lf.\n", gaussTotalMaxDiff);
    printf ("Gauss min difference: %23.16lf.\n", gaussTotalMinDiff);
    printf ("Gauss average difference: %23.16lf.\n\n", gaussTotalAverageDiff);

    printf ("## 4\nAverage gauss elimination time: %dms.\n\n",
            (int) round (totalGauss * 1000.0 / CLOCKS_PER_SEC / RUN_COUNT));

    printf ("## 5\nAverage LUP build time: %dms.\n\n",
            (int) round (totalLUPBuild * 1000.0 / CLOCKS_PER_SEC / RUN_COUNT));

    printf ("## 6\nLUP max difference: %23.16lf.\n", lupTotalMaxDiff);
    printf ("LUP min difference: %23.16lf.\n", lupTotalMinDiff);
    printf ("LUP average difference: %23.16lf.\n\n", lupTotalAverageDiff);

    printf ("## 7\nAverage LUP solve time: %dms.\n\n",
            (int) round (totalLUPSolve * 1000.0 / CLOCKS_PER_SEC / RUN_COUNT));

    printf ("## 8\nCholesky max difference: %23.16lf.\n", choleskyTotalMaxDiff);
    printf ("Cholesky min difference: %23.16lf.\n", choleskyTotalMinDiff);
    printf ("Cholesky average difference: %23.16lf.\n\n", choleskyTotalAverageDiff);

    printf ("## 9\nAverage Cholesky solve time: %dms.\n\n",
            (int) round (totalCholeskySolve * 1000.0 / CLOCKS_PER_SEC / RUN_COUNT));

    printf ("## 10\nRelaxation max difference: %23.16lf.\n", relaxationTotalMaxDiff);
    printf ("Relaxation min difference: %23.16lf.\n", relaxationTotalMinDiff);
    printf ("Relaxation average difference: %23.16lf.\n\n", relaxationTotalAverageDiff);

    printf ("## 11\nAverage relaxation elimination time: %dms.\n\n",
            (int) round (totalRelaxation * 1000.0 / CLOCKS_PER_SEC / RUN_COUNT));

    free (GlobalRand);
    return 0;
}
