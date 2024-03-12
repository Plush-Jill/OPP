#include <math.h>
#include <mpi/mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 15000
#define EPSILON 1E-7
#define TAU 1E-5
#define MAX_ITERATION_COUNT 100000

void initLineCountsAndOffsets(int* lineCountsPerProcess, int* lineOffsets, int processSize);
void initMatrixPart(double* matrixAPart, int lineCountInCurrentProcess, int offsetForCurrentProcess);
void initVectorX(double* vectorX, int size);
void initVectorB(double* vectorB, int size);
double calculateEuclideanVectorNorm(const double* vector, int size);
double calculateSquareRootOfNorm(const double* vector, int size);
void calc_Axb(const double* matrixAPart, const double* vectorX, const double* vectorB, double* vectorAxbPart, int partSize, int partOffset);
void setVectorXtoFofX(const double* vectorAxbPart, const double* vectorX, double* vectorXPart, int partSize, int partOffset);
void printMatrix(const double* matrixAPart, int lineCountInCurrentProcess, int processSize, int rank);
/*void initAllNecessaryArrays(int** lineCounts,
                            int** lineOffsets,
                            double** matrixAPart,
                            double** vectorX,
                            double** vectorB,
                            double** vectorAxbPart,
                            double** vectorXPart,
                            int processCurrentRank,
                            int processTotalSize);*/

int main(int argc, char **argv) {
    int processCurrentRank;
    int processTotalSize;
    int iterationsCount;
    double vectorBEuclideanNorm;
    double accuracy = EPSILON + 1;
    double beginningTime;
    double endingTime;
    int* lineCounts;
    int* lineOffsets;
    double* matrixAPart;
    double* vectorX;
    double* vectorB;
    double* vectorAxbPart;
    double* vectorXPart;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &processTotalSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &processCurrentRank);

    lineCounts = malloc(sizeof(int) * processTotalSize);
    lineOffsets = malloc(sizeof(int) * processTotalSize);
    initLineCountsAndOffsets(lineCounts, lineOffsets, processTotalSize);
    matrixAPart = malloc(sizeof(double) * lineCounts[processCurrentRank] * N);
    vectorX = malloc(sizeof(double) * N);
    vectorB = malloc(sizeof(double) * N);

    initMatrixPart(matrixAPart, lineCounts[processCurrentRank], lineOffsets[processCurrentRank]);
    initVectorX(vectorX, N);
    initVectorB(vectorB, N);
    vectorBEuclideanNorm = 0;
    if (processCurrentRank == 0) {
        printf("All necessary data has been initialized.\nTotal process count: %d\n", processCurrentRank);
        vectorBEuclideanNorm = calculateSquareRootOfNorm(vectorB, N);
    }
    vectorAxbPart = malloc(sizeof(double) * lineCounts[processCurrentRank]);
    vectorXPart = malloc(sizeof(double) * lineCounts[processCurrentRank]);

    beginningTime = MPI_Wtime();
    double vectorAxbEuclideanNorm = 0;
    for (iterationsCount = 0; accuracy > EPSILON && iterationsCount < MAX_ITERATION_COUNT; ++iterationsCount) {

        calc_Axb(matrixAPart, vectorX, vectorB, vectorAxbPart, lineCounts[processCurrentRank], lineOffsets[processCurrentRank]);
        setVectorXtoFofX(vectorAxbPart, vectorX, vectorXPart, lineCounts[processCurrentRank],lineOffsets[processCurrentRank]);

        MPI_Allgatherv(vectorXPart, lineCounts[processCurrentRank], MPI_DOUBLE, vectorX, lineCounts, lineOffsets, MPI_DOUBLE, MPI_COMM_WORLD);

        double vectorAxbPartEuclideanNorm = calculateEuclideanVectorNorm(vectorAxbPart, lineCounts[processCurrentRank]);
        MPI_Reduce(&vectorAxbPartEuclideanNorm, &vectorAxbEuclideanNorm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        if (processCurrentRank == 0) {
            accuracy = sqrt(vectorAxbEuclideanNorm) / vectorBEuclideanNorm;
        }
        MPI_Bcast(&accuracy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    endingTime = MPI_Wtime();

    if (processCurrentRank == 0) {
        if (iterationsCount == MAX_ITERATION_COUNT) {
            printf("Too many iterations\n");
        }else{
            printf("Time: %lf sec\n\n", endingTime - beginningTime);
        }
    }

    free(lineCounts);
    free(lineOffsets);
    free(vectorX);
    free(vectorB);
    free(matrixAPart);
    free(vectorAxbPart);
    free(vectorXPart);

    MPI_Finalize();

    return 0;
}
void printMatrix(const double* matrixAPart, int lineCountInCurrentProcess, int processSize, int rank) {
    for (int processRank = 0; processRank < processSize; ++processRank) {
        if (processRank == rank) {
            for (int i = 0; i < lineCountInCurrentProcess; ++i) {
                for (int j = 0; j < N; ++j) {
                    printf("%lf ", matrixAPart[i * N + j]);
                }
                printf("\n");
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}
void initMatrixPart(double* matrixAPart, int lineCountInCurrentProcess, int offsetForCurrentProcess) {
    for (int i = 0; i < lineCountInCurrentProcess; ++i) {
        for (int j = 0; j < N; ++j) {
            matrixAPart[i * N + j] = 1;
        }
        matrixAPart[i * N + offsetForCurrentProcess + i] = 2;
    }
}
void initVectorX(double* vectorX, int size) {
    memset(vectorX, 0, size * sizeof(double));
}
void initVectorB(double* vectorB, int size) {
    for (int i = 0; i < size; ++i) {
        vectorB[i] = N + 1;
    }
}
void initLineCountsAndOffsets(int* lineCountsPerProcess, int* lineOffsets, int processSize) {
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
void calc_Axb(const double* matrixAPart, const double* vectorX, const double* vectorB, double* vectorAxbPart, int partSize, int partOffset) {
    for (int i = 0; i < partSize; ++i) {
        vectorAxbPart[i] = -vectorB[partOffset + i];
        for (int j = 0; j < N; ++j) {
            vectorAxbPart[i] += matrixAPart[i * N + j] * vectorX[j];
        }
    }
}
void setVectorXtoFofX(const double* vectorAxbPart, const double* vectorX, double* vectorXPart, int partSize, int partOffset) {
    for (int i = 0; i < partSize; ++i) {
        vectorXPart[i] = vectorX[partOffset + i] - TAU * vectorAxbPart[i];
    }
}
double calculateEuclideanVectorNorm(const double* vector, int size) {
    double vectorNorm = 0.0;
    for (int i = 0; i < size; ++i) {
        vectorNorm += vector[i] * vector[i];
    }
    return vectorNorm;
}
double calculateSquareRootOfNorm(const double* vector, int size) {
    return sqrt(calculateEuclideanVectorNorm(vector, size));
}
/*void initAllNecessaryArrays(int** lineCounts,
                            int** lineOffsets,
                            double** matrixAPart,
                            double** vectorX,
                            double** vectorB,
                            double** vectorAxbPart,
                            double** vectorXPart,
                            int processCurrentRank,
                            int processTotalSize){
    *lineCounts = malloc(sizeof(int) * processTotalSize);
    *lineOffsets = malloc(sizeof(int) * processTotalSize);
    initLineCountsAndOffsets(*lineCounts, *lineOffsets, processTotalSize);
    *matrixAPart = malloc(sizeof(double) * (*lineCounts)[processCurrentRank] * N);
    *vectorX = malloc(sizeof(double) * N);
    *vectorB = malloc(sizeof(double) * N);
    initMatrixPart(*matrixAPart, (*lineCounts)[processCurrentRank], (*lineOffsets)[processCurrentRank]);
    initVectorX(*vectorX, N);
    initVectorB(*vectorB, N);
    *vectorAxbPart = malloc(sizeof(double) * (*lineCounts)[processCurrentRank]);
    *vectorXPart = malloc(sizeof(double) * (*lineCounts)[processCurrentRank]);
}*/