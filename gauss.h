#pragma once
#include <stdbool.h>
#include <math.h>

#define FREE_MARKER INFINITY

bool Gauss (double **A, double **B, int *Xi, int rows, int cols, int results);
void PrintGaussSolution (double **sA, double **sB, int *Xi, int rows, int cols, int results);
