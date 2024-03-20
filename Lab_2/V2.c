#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define N 15000
#define EPSILON 1E-7
#define TAU 1E-5
#define MAX_ITERATION_COUNT 100000


void initMatrixA_v2 (double* matrixA){
    for (int i = 0; i < N; ++i){
        for (int j = 0; j < N; ++j){
            matrixA[i * N + j] = 1;
            if (i == j){
                matrixA[i * N + j] = 2;
            }
        }
    }
}
void initVectorX_v2 (double* vectorX){
    memset(vectorX, 0, N * sizeof(double));
}
void initVectorB_v2 (double* vectorB){
    memset(vectorB, N + 1, N * sizeof(double));
}
void calcAxb_v2 (double* vectorAxb, const double* matrixA, const double* vectorX, const double* vectorB, int currentProcessSize){
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            vectorAxb[i] += matrixA[i * N + j] * vectorX[j];
        }
        vectorAxb[i] -= vectorB[i];
    }
}
void setVectorXtoFofX_v2 (double* vectorX, const double* vectorAxb, int sizeInCurrentProcess){
    for (int i = 0; i < N; ++i){
        vectorX[i] -= TAU * (vectorAxb[i]);
    }
}
double calculateEuclideanVectorNorm_v2 (const double* vector, int currentProcessSize){
    double norm = 0;
    for (int i = 0; i < N; ++i){
        norm += vector[i] * vector[i];
    }
    return norm;
}
double calculateSquareRootOfNorm_v2 (const double* vector, int currentProcessSize) {
    return sqrt(calculateEuclideanVectorNorm_v2(vector, currentProcessSize));
}
void printVector_v2 (double* vector){
    printf("(");
    for (int i = 0; i < N; ++i){
        i == N - 1 ? printf("%lf)", vector[i]) : printf("%lf, ", vector[i]);
    }
    //printf(")\n");
}
void initLineCountsAndOffsets_v2(int* lineCountsPerProcess, int* lineOffsets, int processSize) {
    int offset = 0;
    for (int i = 0; i < processSize; ++i) {
        lineCountsPerProcess[i] = N / processSize;

        if (i < N % processSize) {
            ++lineCountsPerProcess[i];
        }
        lineOffsets[i] = offset;
        offset += lineCountsPerProcess[i];
    }
}


int main() {
    int iterationsCount;
    double* vectorX = malloc(sizeof(double) * N);
    double* vectorB = malloc(sizeof(double) * N);
    double* vectorAxb = malloc(sizeof(double) * N);
    double* matrixA = malloc(sizeof(double) * N * N);
    double vectorAxbEuclideanNorm;
    double vectorBEuclideanNorm;
    double accuracy = N;

    int threadID;
    int threadCount = omp_get_max_threads();
    int* lineCounts = calloc(threadCount, sizeof(int));
    int* lineOffsets = calloc(threadCount, sizeof(int));

    initLineCountsAndOffsets_v2(lineCounts, lineOffsets, threadCount);
    initVectorX_v2(vectorX);
    initVectorB_v2(vectorB);
    initMatrixA_v2(matrixA);

    double beginningTime;
    double endingTime;

    beginningTime = omp_get_wtime();
    double vectorAxbNormParts[threadCount];
    vectorBEuclideanNorm = calculateEuclideanVectorNorm_v2(vectorB, N);

    #pragma omp parallel private(thread_id)
    {
        threadID = omp_get_thread_num();

        for (iterationsCount = 0; accuracy >= EPSILON && iterationsCount < MAX_ITERATION_COUNT; ++iterationsCount) {
            calcAxb_v2(vectorAxb + lineOffsets[threadID],
                       matrixA + lineOffsets[threadID] * N,
                       vectorX,
                       vectorB + lineOffsets[threadID],
                       lineCounts[threadID]);

            #pragma omp barrier
            setVectorXtoFofX_v2(vectorX + lineOffsets[threadID], vectorAxb + lineOffsets[threadID],
                                lineCounts[threadID]);

            vectorAxbNormParts[threadID] = calculateEuclideanVectorNorm_v2(vectorAxb + lineOffsets[threadID],
                                                                           lineCounts[threadID]);
            accuracy = 0;

            #pragma omp atomic
            vectorAxbEuclideanNorm += vectorAxbNormParts[threadID];

            #pragma omp barrier
            #pragma omp single
            accuracy = sqrt(vectorAxbEuclideanNorm / vectorBEuclideanNorm);

    }
    }
    endingTime = omp_get_wtime();
    fprintf(stdout, "Time: %lf sec.\n", endingTime - beginningTime);
    printVector_v2(vectorX);


    free(lineOffsets);
    free(lineCounts);
    free(vectorX);
    free(vectorB);
    free(vectorAxb);
    free(matrixA);
}