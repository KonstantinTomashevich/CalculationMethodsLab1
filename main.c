#include <stdio.h>
#include <stdlib.h>

#include "matrixutils.h"
#include "mtwister.h"
#include "gaussjordan.h"

#define abs(a) (((a) > 0.0 ? (a) : (-a)))
#define max(a, b) (((a) > (b) ? (a) : (b)))

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
        double aMax = abs (A[0][0]);
        double antiAMax = abs (antiA[0][0]);

        for (int row = 0; row < MATRIX_SIZE; ++row)
        {
            for (int col = 0; col < MATRIX_SIZE; ++col)
            {
                aMax = max (aMax, abs (A [row] [col]));
                antiAMax = max (antiAMax, abs (antiA [row] [col]));
            }
        }

        printf ("Condition number: %lf.\n", aMax * antiAMax);
    }

    FreeMatrix (copyA, MATRIX_SIZE, MATRIX_SIZE);
    FreeMatrix (antiA, MATRIX_SIZE, MATRIX_SIZE);
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
    FreeMatrix (A, MATRIX_SIZE, MATRIX_SIZE);
    FreeMatrix (X, MATRIX_SIZE, 1);
    FreeMatrix (B, MATRIX_SIZE, 1);
}

int main ()
{
    GlobalRand = malloc (sizeof (MTRand));
    *GlobalRand = SeedRand (1337);

    MainCycle ();
    free (GlobalRand);
    return 0;
}
