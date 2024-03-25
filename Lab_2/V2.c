#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define N 16000
#define EPSILON 1E-10
#define TAU 1E-5
#define MAX_ITERATION_COUNT 1000000

void setLineCountsAndOffsets_v2(int* lineCounts, int* line_offsets, int size, int thread_count);
void initMatrixA_v2(double* matrixA, int size);
void initVectorX_v2(double* vectorX, int size);
void initVectorB_v2(double* vectorB, int size);
double getNormSquare_v2(const double* vector, int size);
void calcVectorAXB_v2(const double* matrixA, const double* vectorX, const double* vectorB, double* vectorAXB, int size);
void calcNextVectorX_v2(const double* vectorAXB, double* vectorX, double tau, int size);
bool checkResultVector_v200 (const double* vectorX) {
    double preEpsilon = 1 - EPSILON;
    double postEpsilon = 1 + EPSILON;

    for (int i = 0; i < N; ++i) {
        if (!(vectorX[i] < postEpsilon && vectorX[i] > preEpsilon)) {
            return false;
        }
    }

    return true;
}

int main() {
    double* matrixA = malloc(sizeof(double) * N * N);
    double* vectorX = malloc(sizeof(double) * N);
    double* vectorB = malloc(sizeof(double) * N);
    double* vectorAXB = malloc(sizeof(double) * N);

    initMatrixA_v2(matrixA, N);
    initVectorX_v2(vectorX, N);
    initVectorB_v2(vectorB, N);


    double vectorBNorm = sqrt(getNormSquare_v2(vectorB, N));


    int i = 0;
    double accuracy = EPSILON;
    int currentThreadCount;
    int* lineCounts;
    int* lineOffsets;

    double beginningTime = omp_get_wtime();

    #pragma omp parallel shared(lineCounts, lineOffsets)
    {
        currentThreadCount = omp_get_num_threads();
        int threadID = omp_get_thread_num();

        lineCounts = malloc(sizeof(int) * currentThreadCount);
        lineOffsets = malloc(sizeof(int) * currentThreadCount);
        setLineCountsAndOffsets_v2(lineCounts, lineOffsets, N, currentThreadCount);

        #pragma omp single
        {
            fprintf(stdout, "currentThreadCount: %d\n", currentThreadCount);
        };
        for (; accuracy >= EPSILON && i < MAX_ITERATION_COUNT; ++i) {
            calcVectorAXB_v2(matrixA + lineOffsets[threadID] * N, vectorX, vectorB + lineOffsets[threadID],
                             vectorAXB + lineOffsets[threadID], lineCounts[threadID]);
            #pragma omp barrier

            calcNextVectorX_v2(vectorAXB + lineOffsets[threadID], vectorX + lineOffsets[threadID], TAU,
                               lineCounts[threadID]);

            #pragma omp single
            accuracy = 0;

            #pragma omp atomic
            accuracy += getNormSquare_v2(vectorAXB + lineOffsets[threadID], lineCounts[threadID]);
            #pragma omp barrier

            #pragma omp single
            accuracy = sqrt(accuracy) / vectorBNorm;
        }
    }

    double endingTime = omp_get_wtime();

    fprintf(stdout, "Time = %.10lf sec\n", endingTime - beginningTime);
    if (checkResultVector_v200(vectorX)) {
        printf("Vector is calculated correctly.\n");
    } else {
        printf("Vector is calculated incorrectly.\n");
    }

    free(lineCounts);
    free(lineOffsets);
    free(matrixA);
    free(vectorX);
    free(vectorB);
    free(vectorAXB);

    return 0;
}

void setLineCountsAndOffsets_v2(int* lineCounts, int* line_offsets, int size, int thread_count) {
    int offset = 0;
    for (int i = 0; i < thread_count; ++i)
    {
        lineCounts[i] = size / thread_count;

        if (i < size % thread_count) {
            ++lineCounts[i];
        }
        line_offsets[i] = offset;
        offset += lineCounts[i];
    }
}
void initMatrixA_v2(double* matrixA, int size) {
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            matrixA[i * size + j] = 1;
        }
        matrixA[i * size + i] = 2;
    }
}
void initVectorX_v2(double* vectorX, int size) {
    for (int i = 0; i < size; ++i) {
        vectorX[i] = 0;
    }
}
void initVectorB_v2(double* vectorB, int size) {
    for (int i = 0; i < size; ++i)
        vectorB[i] = N + 1;
}
double getNormSquare_v2(const double* vector, int size) {
    double res = 0.0;
    for (int i = 0; i < size; ++i)
        res += vector[i] * vector[i];

    return res;
}
void calcVectorAXB_v2(const double* matrixA, const double* vectorX, const double* vectorB, double* vectorAXB, int size) {
    for (int i = 0; i < size; ++i) {
        vectorAXB[i] = -vectorB[i];
        for (int j = 0; j < N; ++j)
            vectorAXB[i] += matrixA[i * N + j] * vectorX[j];
    }
}
void calcNextVectorX_v2(const double* vectorAXB, double* vectorX, double tau, int size) {
    for (int i = 0; i < size; ++i) {
        vectorX[i] -= tau * vectorAXB[i];
    }
}
