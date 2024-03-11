#include <math.h>
#include <mpi/mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 10000
#define EPSILON 1E-7
#define TAU 1E-5
#define MAX_ITERATION_COUNT 100000

void setMatrixPart(int* lineCountsPerProcess, int* lineOffsets, int processSize);
void initMatrixChunks(double* matrixAPart, int lineCountInCurrentProcess, int offsetForCurrentProcess);
void initVectorX(double* vectorX, int size);
void initVectorB(double* vectorB, int size);
double calculateEuclideanVectorNorm(const double* vector, int size);
double calculateSquareRootOfNorm(const double* vector, int size);
void calc_Axb(const double* matrixAPart, const double* vectorX, const double* vectorB, double* matrixAxbPart, int chunkSize, int chunkOffset);
void calc_FofX(const double* matrixAxbPart, const double* vectorX, double* vectorXPart, int chunkSize, int chunkOffset);
void printMatrix(const double* matrixAPart, int lineCountInCurrentProcess, int processSize, int rank);

int main(int argc, char **argv) {
    int rank;
    int size;
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
    double* AxbPart;
    double* vectorXPart;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    lineCounts = malloc(sizeof(int) * size);
    lineOffsets = malloc(sizeof(int) * size);
    setMatrixPart(lineCounts, lineOffsets, size);
    matrixAPart = malloc(sizeof(double) * lineCounts[rank] * N);
    vectorX = malloc(sizeof(double) * N);
    vectorB = malloc(sizeof(double) * N);

    initMatrixChunks(matrixAPart, lineCounts[rank], lineOffsets[rank]);
    initVectorX(vectorX, N);
    initVectorB(vectorB, N);
    vectorBEuclideanNorm = 0;
    if (rank == 0) {
        printf("START\n");
        vectorBEuclideanNorm = calculateSquareRootOfNorm(vectorB, N);
    }
    AxbPart = malloc(sizeof(double) * lineCounts[rank]);
    vectorXPart = malloc(sizeof(double) * lineCounts[rank]);

    beginningTime = MPI_Wtime();
    double AxbEuclideanNorm = 0;
    for (iterationsCount = 0; accuracy > EPSILON && iterationsCount < MAX_ITERATION_COUNT; ++iterationsCount) {

        calc_Axb(matrixAPart, vectorX, vectorB, AxbPart, lineCounts[rank], lineOffsets[rank]);
        calc_FofX(AxbPart, vectorX, vectorXPart, lineCounts[rank], lineOffsets[rank]);

        MPI_Allgatherv(vectorXPart, lineCounts[rank], MPI_DOUBLE, vectorX, lineCounts, lineOffsets, MPI_DOUBLE, MPI_COMM_WORLD);

        double AxbPartEuclideanNorm = calculateEuclideanVectorNorm(AxbPart, lineCounts[rank]);
        MPI_Reduce(&AxbPartEuclideanNorm, &AxbEuclideanNorm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            accuracy = sqrt(AxbEuclideanNorm) / vectorBEuclideanNorm;
        }
        MPI_Bcast(&accuracy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    endingTime = MPI_Wtime();

    if (rank == 0) {
        if (iterationsCount == MAX_ITERATION_COUNT) {
            printf("Too many iterations\n");
        }else{
            printf("Time: %lf sec\n\n", endingTime - beginningTime);
        }
    }
    if (N <= 15){
        printMatrix(matrixAPart, lineCounts[rank], size, rank);
    }

    free(lineCounts);
    free(lineOffsets);
    free(vectorX);
    free(vectorB);
    free(matrixAPart);
    free(AxbPart);
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
void initMatrixChunks(double* matrixAPart, int lineCountInCurrentProcess, int offsetForCurrentProcess) {
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
void setMatrixPart(int* lineCountsPerProcess, int* lineOffsets, int processSize) {
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
void calc_Axb(const double* matrixAPart, const double* vectorX, const double* vectorB, double* matrixAxbPart, int chunkSize, int chunkOffset) {
    for (int i = 0; i < chunkSize; ++i) {
        matrixAxbPart[i] = -vectorB[chunkOffset + i];
        for (int j = 0; j < N; ++j) {
            matrixAxbPart[i] += matrixAPart[i * N + j] * vectorX[j];
        }
    }
}
void calc_FofX(const double* matrixAxbPart, const double* vectorX, double* vectorXPart, int chunkSize, int chunkOffset) {
    for (int i = 0; i < chunkSize; ++i) {
        vectorXPart[i] = vectorX[chunkOffset + i] - TAU * matrixAxbPart[i];
    }
}
double calculateEuclideanVectorNorm(const double* vector, int size) {
    double vectorNorm = 0.0;
    for (int i = 0; i < size; ++i) {
        vectorNorm += vector[i] * vector[i];
    }
    return vectorNorm;
}
double calculateSquareRootOfNorm(const double* vector, int size){
    return sqrt(calculateEuclideanVectorNorm(vector, size));
}
