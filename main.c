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
#include "householder.h"
#include "minquads.h"
#include "gmres.h"
#include "gmresarnoldi.h"

#define RUN_COUNT 100
extern MTRand *GlobalRand;
double minConditionNumber = INFINITY;
double maxConditionNumber = 0.0;
double averageMatrixElement = 0.0;

clock_t totalGaussJordan = 0;

double gaussTotalMaxNormal = 0.0;
double gaussTotalMinNormal = INFINITY;
double gaussTotalAverageNormal = 0.0;

clock_t totalGauss = 0;

clock_t totalLUPBuild = 0;

double lupTotalMaxNormal = 0.0;
double lupTotalMinNormal = INFINITY;
double lupTotalAverageNormal = 0.0;

clock_t totalLUPSolve = 0;

double choleskyTotalMaxNormal = 0.0;
double choleskyTotalMinNormal = INFINITY;
double choleskyTotalAverageNormal = 0.0;

clock_t totalCholeskySolve = 0;

double relaxationTotalMaxNormal = 0.0;
double relaxationTotalMinNormal = INFINITY;
double relaxationTotalAverageNormal = 0.0;

clock_t totalRelaxation = 0;

double householderTotalMaxNormal = 0.0;
double householderTotalMinNormal = INFINITY;
double householderTotalAverageNormal = 0.0;

clock_t totalHouseholder = 0;

double minquadsTotalMaxNormal = 0.0;
double minquadsTotalMinNormal = INFINITY;
double minquadsTotalAverageNormal = 0.0;

clock_t totalMinQuads = 0;

double gmresTotalMaxNormal = 0.0;
double gmresTotalMinNormal = INFINITY;
double gmresTotalAverageNormal = 0.0;

clock_t totalGMRES = 0;

double gmresArnoldiTotalMaxNormal = 0.0;
double gmresArnoldiTotalMinNormal = INFINITY;
double gmresArnoldiTotalAverageNormal = 0.0;

clock_t totalGMRESArnoldi = 0;

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
        double aMax = 0.0;
        double antiAMax = 0.0;
        double average = 0.0;

        for (int row = 0; row < MATRIX_SIZE; ++row)
        {
            double aSum = 0.0;
            double antiASum = 0.0;

            for (int col = 0; col < MATRIX_SIZE; ++col)
            {
                average += A[row][col] / (MATRIX_SIZE * MATRIX_SIZE);
                aSum += fabs (A[row][col]);
                antiASum += fabs (antiA[row][col]);
            }

            aMax = fmax (aMax, aSum);
            antiAMax = fmax (antiAMax, antiASum);
        }

        minConditionNumber = fmin (minConditionNumber, aMax * antiAMax);
        maxConditionNumber = fmax (maxConditionNumber, aMax * antiAMax);
        averageMatrixElement += average / RUN_COUNT;
        printf ("Condition number: %lf.\n", aMax * antiAMax);
    }

    FreeMatrix (copyA, MATRIX_SIZE, MATRIX_SIZE);
    FreeMatrix (antiA, MATRIX_SIZE, MATRIX_SIZE);
}

void FindGaussSolutionAndPrintNormal (double **A, double **B, double **X)
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
        double normal = 0.0;

        for (int index = 0; index < MATRIX_SIZE; ++index)
        {
            double currentNormal = fabs (X[Xi[index]][0] - copyB[index][0]);
            normal += currentNormal * currentNormal;
        }

        normal = sqrt (normal);
        printf ("Gauss normal (quadric): %23.16lf.\n", normal);

        gaussTotalMaxNormal = fmax (gaussTotalMaxNormal, normal);
        gaussTotalMinNormal = fmin (gaussTotalMinNormal, normal);
        gaussTotalAverageNormal += normal / RUN_COUNT;
    }

    FreeMatrix (copyA, MATRIX_SIZE, MATRIX_SIZE);
    FreeMatrix (copyB, MATRIX_SIZE, 1);
    free (Xi);
}

void FindLUPSolutionAndPrintNormal (double **A, double **B, double **X)
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
            double normal = 0.0;

            for (int index = 0; index < MATRIX_SIZE; ++index)
            {
                double currentNormal = fabs (X[index][0] - builtX[index][0]);
                normal += currentNormal * currentNormal;
            }

            normal = sqrt (normal);
            printf ("LUP normal (quadric): %23.16lf.\n", normal);

            lupTotalMaxNormal = fmax (lupTotalMaxNormal, normal);
            lupTotalMinNormal = fmin (lupTotalMinNormal, normal);
            lupTotalAverageNormal += normal / RUN_COUNT;
            FreeMatrix (builtX, MATRIX_SIZE, MATRIX_SIZE);
        }
    }

    FreeMatrix (LU, MATRIX_SIZE, MATRIX_SIZE);
    FreeMatrix (copyB, MATRIX_SIZE, 1);
}

void FindCholeskySolutionAndPrintNormal (double **A, double **B, double **X)
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
            double normal = 0.0;

            for (int index = 0; index < MATRIX_SIZE; ++index)
            {
                double currentNormal = fabs (X[index][0] - builtX[index][0]);
                normal += currentNormal * currentNormal;
            }

            normal = sqrt (normal);
            printf ("Cholesky normal (quadric): %23.16lf.\n", normal);

            choleskyTotalMaxNormal = fmax (choleskyTotalMaxNormal, normal);
            choleskyTotalMinNormal = fmin (choleskyTotalMinNormal, normal);
            choleskyTotalAverageNormal += normal / RUN_COUNT;
            FreeMatrix (builtX, MATRIX_SIZE, MATRIX_SIZE);
        }
    }

    FreeMatrix (LT, MATRIX_SIZE, MATRIX_SIZE);
    FreeMatrix (copyB, MATRIX_SIZE, 1);
}

void FindRelaxationSolutionAndPrintNormal (double **A, double **B, double **X)
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
        double normal = 0.0;

        for (int index = 0; index < MATRIX_SIZE; ++index)
        {
            double currentNormal = fabs (builtX[index][0] - X[index][0]);
            normal += currentNormal * currentNormal;
        }

        normal = sqrt (normal);
        printf ("Relaxation normal (quadric): %23.16lf.\n", normal);

        relaxationTotalMaxNormal = fmax (relaxationTotalMaxNormal, normal);
        relaxationTotalMinNormal = fmin (relaxationTotalMinNormal, normal);
        relaxationTotalAverageNormal += normal / RUN_COUNT;
    }

    FreeMatrix (copyA, MATRIX_SIZE, MATRIX_SIZE);
    FreeMatrix (copyB, MATRIX_SIZE, 1);
}

void FindHouseholderSolutionAndPrintNormal (double **A, double **B, double **X)
{
    double **copyA = CopyMatrix (A, MATRIX_SIZE, MATRIX_SIZE);
    double **copyB = CopyMatrix (B, MATRIX_SIZE, 1);
    clock_t begin = clock ();

    if (!SolveHouseholder (copyA, MATRIX_SIZE, copyB, 1))
    {
        printf ("Unable to solve system!\n");
    }
    else
    {
        totalHouseholder += clock () - begin;
        double normal = 0.0;

        for (int index = 0; index < MATRIX_SIZE; ++index)
        {
            double currentNormal = fabs (copyB[index][0] - X[index][0]);
            normal += currentNormal * currentNormal;
        }

        normal = sqrt (normal);
        printf ("Householder normal (quadric): %23.16lf.\n", normal);

        householderTotalMaxNormal = fmax (householderTotalMaxNormal, normal);
        householderTotalMinNormal = fmin (householderTotalMinNormal, normal);
        householderTotalAverageNormal += normal / RUN_COUNT;
    }

    FreeMatrix (copyA, MATRIX_SIZE, MATRIX_SIZE);
    FreeMatrix (copyB, MATRIX_SIZE, 1);
}

void FindMinQuadsSolutionAndPrintNormal (double **A, double **B, double **X)
{
    double **copyA = CopyMatrix (A, MATRIX_SIZE, MATRIX_SIZE);
    double **copyB = CopyMatrix (B, MATRIX_SIZE, 1);
    double **builtX;
    int colsCount = N * 20;
    clock_t begin = clock ();

    if (!SolveMinQuads (copyA, MATRIX_SIZE, colsCount, copyB, 1, &builtX))
    {
        printf ("Unable to solve system!\n");
    }
    else
    {
        totalMinQuads += clock () - begin;
        double normal = 0.0;

        for (int index = 0; index < colsCount; ++index)
        {
            double currentNormal = fabs (builtX[index][0] - X[index][0]);
            normal += currentNormal * currentNormal;
        }

        normal = sqrt (normal);
        printf ("MinQuads normal (quadric): %23.16lf.\n", normal);

        minquadsTotalMaxNormal = fmax (minquadsTotalMaxNormal, normal);
        minquadsTotalMinNormal = fmin (minquadsTotalMinNormal, normal);
        minquadsTotalAverageNormal += normal / RUN_COUNT;
    }

    FreeMatrix (copyA, MATRIX_SIZE, MATRIX_SIZE);
    FreeMatrix (copyB, MATRIX_SIZE, 1);
}

void FindGMRESSolutionAndPrintNormal (double **A, double **B, double **X)
{
    double **copyA = CopyMatrix (A, MATRIX_SIZE, MATRIX_SIZE);
    double **copyB = CopyMatrix (B, MATRIX_SIZE, 1);
    double **builtX;
    clock_t begin = clock ();

    if (!SolveGMRES (copyA, MATRIX_SIZE, copyB, &builtX))
    {
        printf ("Unable to solve system!\n");
    }
    else
    {
        totalGMRES += clock () - begin;
        double normal = 0.0;

        for (int index = 0; index < MATRIX_SIZE; ++index)
        {
            double currentNormal = fabs (builtX[index][0] - X[index][0]);
            normal += currentNormal * currentNormal;
        }

        normal = sqrt (normal);
        printf ("GMRES normal (quadric): %23.16lf.\n", normal);

        gmresTotalMaxNormal = fmax (gmresTotalMaxNormal, normal);
        gmresTotalMinNormal = fmin (gmresTotalMinNormal, normal);
        gmresTotalAverageNormal += normal / RUN_COUNT;
    }

    FreeMatrix (copyA, MATRIX_SIZE, MATRIX_SIZE);
    FreeMatrix (copyB, MATRIX_SIZE, 1);
}

void FindGMRESArnoldiSolutionAndPrintNormal (double **A, double **B, double **X)
{
    double **copyA = CopyMatrix (A, MATRIX_SIZE, MATRIX_SIZE);
    double **copyB = CopyMatrix (B, MATRIX_SIZE, 1);
    double **builtX;
    clock_t begin = clock ();

    if (!SolveGMRESArnoldi (copyA, MATRIX_SIZE, copyB, &builtX))
    {
        printf ("Unable to solve system!\n");
    }
    else
    {
        totalGMRESArnoldi += clock () - begin;
        double normal = 0.0;

        for (int index = 0; index < MATRIX_SIZE; ++index)
        {
            double currentNormal = fabs (builtX[index][0] - X[index][0]);
            normal += currentNormal * currentNormal;
        }

        normal = sqrt (normal);
        printf ("GMRESArnoldi normal (quadric): %23.16lf.\n", normal);

        gmresArnoldiTotalMaxNormal = fmax (gmresArnoldiTotalMaxNormal, normal);
        gmresArnoldiTotalMinNormal = fmin (gmresArnoldiTotalMinNormal, normal);
        gmresArnoldiTotalAverageNormal += normal / RUN_COUNT;
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
    FindGaussSolutionAndPrintNormal (A, B, X);
    FindLUPSolutionAndPrintNormal (A, B, X);
    FindCholeskySolutionAndPrintNormal (A, B, X);
    FindRelaxationSolutionAndPrintNormal (A, B, X);
    FindHouseholderSolutionAndPrintNormal (A, B, X);
    FindMinQuadsSolutionAndPrintNormal (A, B, X);
    FindGMRESSolutionAndPrintNormal (A, B, X);
    FindGMRESArnoldiSolutionAndPrintNormal (A, B, X);

    FreeMatrix (A, MATRIX_SIZE, MATRIX_SIZE);
    FreeMatrix (X, MATRIX_SIZE, 1);
    FreeMatrix (B, MATRIX_SIZE, 1);
}

void PrintReport (FILE *output)
{
    fprintf (output, "\n### Report ###\n");
    fprintf (output, "## 1\nMax condition number: %23.16lf.\n", maxConditionNumber);
    fprintf (output, "Min condition number: %23.16lf.\n", minConditionNumber);
    fprintf (output, "Average matrix element: %23.16lf.\n\n", averageMatrixElement);

    fprintf (output, "## 2\nAverage A^-1 calculation time: %dms.\n\n",
            (int) round (totalGaussJordan * 1000.0 / CLOCKS_PER_SEC / RUN_COUNT));

    fprintf (output, "## 3\nGauss max normal (quadric): %23.16lf.\n", gaussTotalMaxNormal);
    fprintf (output, "Gauss min normal (quadric): %23.16lf.\n", gaussTotalMinNormal);
    fprintf (output, "Gauss average normal (quadric): %23.16lf.\n\n", gaussTotalAverageNormal);

    fprintf (output, "## 4\nAverage gauss elimination time: %dms.\n\n",
            (int) round (totalGauss * 1000.0 / CLOCKS_PER_SEC / RUN_COUNT));

    fprintf (output, "## 5\nAverage LUP build time: %dms.\n\n",
            (int) round (totalLUPBuild * 1000.0 / CLOCKS_PER_SEC / RUN_COUNT));

    fprintf (output, "## 6\nLUP max normal (quadric): %23.16lf.\n", lupTotalMaxNormal);
    fprintf (output, "LUP min normal (quadric): %23.16lf.\n", lupTotalMinNormal);
    fprintf (output, "LUP average normal (quadric): %23.16lf.\n\n", lupTotalAverageNormal);

    fprintf (output, "## 7\nAverage LUP solve time: %dms.\n\n",
            (int) round (totalLUPSolve * 1000.0 / CLOCKS_PER_SEC / RUN_COUNT));

    fprintf (output, "## 8\nCholesky max normal (quadric): %23.16lf.\n", choleskyTotalMaxNormal);
    fprintf (output, "Cholesky min normal (quadric): %23.16lf.\n", choleskyTotalMinNormal);
    fprintf (output, "Cholesky average normal (quadric): %23.16lf.\n\n", choleskyTotalAverageNormal);

    fprintf (output, "## 9\nAverage Cholesky solve time: %dms.\n\n",
            (int) round (totalCholeskySolve * 1000.0 / CLOCKS_PER_SEC / RUN_COUNT));

    fprintf (output, "## 10\nRelaxation max normal (quadric): %23.16lf.\n", relaxationTotalMaxNormal);
    fprintf (output, "Relaxation min normal (quadric): %23.16lf.\n", relaxationTotalMinNormal);
    fprintf (output, "Relaxation average normal (quadric): %23.16lf.\n\n", relaxationTotalAverageNormal);

    fprintf (output, "## 11\nAverage relaxation elimination time: %dms.\n\n",
            (int) round (totalRelaxation * 1000.0 / CLOCKS_PER_SEC / RUN_COUNT));

    fprintf (output, "## 12\nHouseholder max normal (quadric): %23.16lf.\n", householderTotalMaxNormal);
    fprintf (output, "Householder min normal (quadric): %23.16lf.\n", householderTotalMinNormal);
    fprintf (output, "Householder average normal (quadric): %23.16lf.\n\n", householderTotalAverageNormal);

    fprintf (output, "## 13\nAverage householder elimination time: %dms.\n\n",
            (int) round (totalHouseholder * 1000.0 / CLOCKS_PER_SEC / RUN_COUNT));

    fprintf (output, "## 14\nMinQuads max normal (quadric): %23.16lf.\n", minquadsTotalMaxNormal);
    fprintf (output, "MinQuads min normal (quadric): %23.16lf.\n", minquadsTotalMinNormal);
    fprintf (output, "MinQuads average normal (quadric): %23.16lf.\n\n", minquadsTotalAverageNormal);

    fprintf (output, "## 15\nAverage minquads elimination time: %dms.\n\n",
            (int) round (totalMinQuads * 1000.0 / CLOCKS_PER_SEC / RUN_COUNT));

    fprintf (output, "## 16\nGMRES max normal (quadric): %23.16lf.\n", gmresTotalMaxNormal);
    fprintf (output, "GMRES min normal (quadric): %23.16lf.\n", gmresTotalMinNormal);
    fprintf (output, "GMRES average normal (quadric): %23.16lf.\n\n", gmresTotalAverageNormal);

    fprintf (output, "## 17\nAverage gmres elimination time: %dms.\n\n",
            (int) round (totalGMRES * 1000.0 / CLOCKS_PER_SEC / RUN_COUNT));

    fprintf (output, "## 18\nGMRESArnoldi max normal (quadric): %23.16lf.\n", gmresArnoldiTotalMaxNormal);
    fprintf (output, "GMRESArnoldi min normal (quadric): %23.16lf.\n", gmresArnoldiTotalMinNormal);
    fprintf (output, "GMRESArnoldi average normal (quadric): %23.16lf.\n\n", gmresArnoldiTotalAverageNormal);

    fprintf (output, "## 19\nAverage gmres arnoldi elimination time: %dms.\n\n",
            (int) round (totalGMRESArnoldi * 1000.0 / CLOCKS_PER_SEC / RUN_COUNT));
}

int main ()
{
    // TODO: After finishing all algos, speed up them by addition of "start from row/col" parameter to util functions.
    GlobalRand = malloc (sizeof (MTRand));
    *GlobalRand = SeedRand (1492);

    for (int index = 0; index < RUN_COUNT; ++index)
    {
        printf ("\n### Run %d ### \n", index);
        MainCycle ();
    }

    PrintReport (stdout);
    FILE *report = fopen ("report.txt", "w");
    PrintReport (report);
    fclose (report);

    free (GlobalRand);
    return 0;
}
