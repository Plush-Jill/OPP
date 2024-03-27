#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>


#define N 16000
#define EPSILON 1E-10
#define TAU 1E-5
#define MAX_ITERATION_COUNT 1000000

void initMatrixA_v1(double* matrixA, int size);
void initVectorX_v1(double* vectorX, int size);
void initVectorB_v1(double* vectorB, int size);
double getNormSquare_v1(const double* vector, int size);
void calcVectorAXB_v1(const double* matrixA, const double* vectorX, const double* vectorB, double* vectorAXB, int size);
void calcNextVectorX_v1(const double* vectorAXB, double* vectorX, double tau, int size);
bool checkResultVector_v1 (const double* vectorX);

int main() {
    double* matrixA = malloc(sizeof(double) * N * N);
    double* vectorX = malloc(sizeof(double) * N);
    double* vectorB = malloc(sizeof(double) * N);
    double* vectorAXB = malloc(sizeof(double) * N);

    initMatrixA_v1(matrixA, N);
    initVectorX_v1(vectorX, N);
    initVectorB_v1(vectorB, N);

    double vectorBNorm = sqrt(getNormSquare_v1(vectorB, N));

    fprintf(stdout, "Threads count: %d\n", omp_get_max_threads());

    double beginningTime = omp_get_wtime();

    int i;
    double accuracy = EPSILON + 1;
    for (i = 0; accuracy > EPSILON && i < MAX_ITERATION_COUNT; ++i) {
        calcVectorAXB_v1(matrixA, vectorX, vectorB, vectorAXB, N);
        calcNextVectorX_v1(vectorAXB, vectorX, TAU, N);
        accuracy = sqrt(getNormSquare_v1(vectorAXB, N)) / vectorBNorm;
    }

    double endingTime = omp_get_wtime();


    fprintf(stdout, "Time = %.10lf sec\n", endingTime - beginningTime);
    if (checkResultVector_v1(vectorX)) {
        printf("Vector is calculated correctly.\n");
    } else {
        printf("Vector is calculated incorrectly.\n");
    }

    free(matrixA);
    free(vectorX);
    free(vectorB);
    free(vectorAXB);

    return 0;
}

void initMatrixA_v1(double* matrixA, int size) {
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            matrixA[i * size + j] = 1;
        }
        matrixA[i * size + i] = 2;
    }
}
void initVectorX_v1(double* vectorX, int size) {
    for (int i = 0; i < size; ++i) {
        vectorX[i] = 0;
    }
}
void initVectorB_v1(double* vectorB, int size) {
    for (int i = 0; i < size; ++i) {
        vectorB[i] = N + 1;
    }
}
double getNormSquare_v1(const double* vector, int size) {
    double normSquare = 0;

    #pragma omp parallel for schedule(runtime) reduction(+:normSquare)
    for (int i = 0; i < size; ++i) {
        normSquare += vector[i] * vector[i];
    }
    return normSquare;
}
void calcVectorAXB_v1(const double* matrixA, const double* vectorX, const double* vectorB, double* vectorAXB, int size) {
    #pragma omp parallel for schedule(runtime)
    for (int i = 0; i < size; ++i) {
        vectorAXB[i] = -vectorB[i];
        for (int j = 0; j < N; ++j)
            vectorAXB[i] += matrixA[i * N + j] * vectorX[j];
    }
}
void calcNextVectorX_v1(const double* vectorAXB, double* vectorX, double tau, int size) {
    #pragma omp parallel for schedule(runtime)
    for (int i = 0; i < size; ++i) {
        vectorX[i] -= tau * vectorAXB[i];
    }


}
bool checkResultVector_v1 (const double* vectorX) {
    double preEpsilon = 1 - EPSILON;
    double postEpsilon = 1 + EPSILON;

    for (int i = 0; i < N; ++i) {
        if (!(vectorX[i] < postEpsilon && vectorX[i] > preEpsilon)) {
            return false;
        }
    }

    return true;
}