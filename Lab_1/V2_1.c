#include <math.h>
#include <mpi/mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define N 10000
#define EPSILON 1E-7
#define TAU 1E-5
#define MAX_ITERATION_COUNT 100000

void setMatrixPart(int* lineCountsPerProcess, int* lineOffsets, int processCount);
void initMatrixAParts(double* matrixAPart, int lineCountInCurrentProcess, int processRank);
void initVectorXParts(double* vectorXPart, int size);
void initVectorBParts(double* vectorBPart, int size);
double calculateEuclideanVectorNorm(const double* vector, int size);
void calculate_Axb(const double* matrixAPart, const double* vectorXPart, const double* vectorBPart, double* replaceVectorXPart,
                   double* matrixAxbPart, int* lineCounts, const int* lineOffsets, int processRank, int processCount);
void calculate_FOfX(const double* matrixAxb, double* vectorXPart, int partSize);
void copyVector(double* dest, const double* src, int size);

int main(int argc, char** argv) {
    int processCurrentRank;
    int processTotalSize;
    int iterationsCount;
    double vectorBPartNorm;
    double vectorBNorm;
    double vectorXPartNorm;
    double vectorXNorm;
    double vectorAxbNorm;
    double accuracy = EPSILON + 1;
    double beginningTime;
    double endingTime;
    int* lineCounts;
    int* lineOffsets;
    double* vectorXPart;
    double* vectorBPart;
    double* matrixAPart;
    double* vectorAxbPart;
    double* replaceVectorXParts;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &processTotalSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &processCurrentRank);

    lineCounts = malloc(sizeof(int) * processTotalSize);
    lineOffsets = malloc(sizeof(int) * processTotalSize);
    setMatrixPart(lineCounts, lineOffsets, processTotalSize);

    vectorXPart = malloc(sizeof(double) * lineCounts[processCurrentRank]);
    vectorBPart = malloc(sizeof(double) * lineCounts[processCurrentRank]);
    matrixAPart = malloc(sizeof(double) * lineCounts[processCurrentRank] * N);
    initVectorXParts(vectorXPart, lineCounts[processCurrentRank]);
    initVectorBParts(vectorBPart, lineCounts[processCurrentRank]);
    initMatrixAParts(matrixAPart, lineCounts[processCurrentRank], lineOffsets[processCurrentRank]);

    vectorBPartNorm = calculateEuclideanVectorNorm(vectorBPart, lineCounts[processCurrentRank]);
    MPI_Reduce(&vectorBPartNorm, &vectorBNorm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (processCurrentRank == 0) {
        vectorBNorm = sqrt(vectorBNorm);
    }
    vectorAxbPart = malloc(sizeof(double) * lineCounts[processCurrentRank]);
    replaceVectorXParts = malloc(sizeof(double) * lineCounts[0]);


    beginningTime = MPI_Wtime();
    for (iterationsCount = 0; accuracy > EPSILON && iterationsCount < MAX_ITERATION_COUNT; ++iterationsCount) {
        calculate_Axb(matrixAPart, vectorXPart, vectorBPart, replaceVectorXParts, vectorAxbPart, lineCounts, lineOffsets, processCurrentRank,
                      processTotalSize);

        calculate_FOfX(vectorAxbPart, vectorXPart, lineCounts[processCurrentRank]);
        vectorAxbNorm = calculateEuclideanVectorNorm(vectorAxbPart, lineCounts[processCurrentRank]);

        MPI_Reduce(&vectorAxbNorm, &accuracy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (processCurrentRank == 0) {
            accuracy = sqrt(accuracy) / vectorBNorm;
        }
        MPI_Bcast(&accuracy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    endingTime = MPI_Wtime();

    vectorXPartNorm = calculateEuclideanVectorNorm(vectorXPart, lineCounts[processCurrentRank]);
    MPI_Reduce(&vectorXPartNorm, &vectorXNorm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (processCurrentRank == 0) {
        if (iterationsCount == MAX_ITERATION_COUNT) {
            fprintf(stderr, "Too many iterations\n");
        }else{
            printf("Norm: %lf\n", sqrt(vectorXNorm));
            printf("Time: %lf sec\n", endingTime - beginningTime);
        }
    }

    free(lineCounts);
    free(lineOffsets);
    free(vectorXPart);
    free(vectorBPart);
    free(matrixAPart);
    free(vectorAxbPart);
    MPI_Finalize();

    return 0;
}

void copyVector(double* dest, const double* src, int size) {
    for (int i = 0; i < size; ++i) {
        dest[i] = src[i];
    }
}
void initMatrixAParts(double* matrixAPart, int lineCountInCurrentProcess, int processRank) {
    for (int i = 0; i < lineCountInCurrentProcess; ++i) {
        for (int j = 0; j < N; ++j) {
            matrixAPart[i * N + j] = 1;
        }
        matrixAPart[i * N + processRank + i] = 2;
    }
}
void initVectorXParts(double* vectorXPart, int size) {
    for (int i = 0; i < size; ++i) {
        vectorXPart[i] = 0;
    }
}
void initVectorBParts(double* vectorBPart, int size) {
    for (int i = 0; i < size; ++i) {
        vectorBPart[i] = N + 1;
    }
}
void setMatrixPart(int* lineCountsPerProcess, int* lineOffsets, int processCount) {
    int offset = 0;
    for (int i = 0; i < processCount; ++i) {
        lineCountsPerProcess[i] = N / processCount;

        if (i < N % processCount) {
            ++lineCountsPerProcess[i];
        }
        lineOffsets[i] = offset;
        offset += lineCountsPerProcess[i];
    }
}
double calculateEuclideanVectorNorm(const double* vector, int size) {
    double norm = 0.0;
    for (int i = 0; i < size; ++i) {
        norm += vector[i] * vector[i];
    }
    return norm;
}
void calculate_FOfX(const double* matrixAxb, double* vectorXPart, int partSize) {
    for (int i = 0; i < partSize; ++i) {
        vectorXPart[i] -= TAU * matrixAxb[i];
    }
}
void calculate_Axb(const double* matrixAPart, const double* vectorXPart, const double* vectorBPart, double* replaceVectorXPart,
                   double* matrixAxbPart, int* lineCounts, const int* lineOffsets, int processRank, int processCount) {
    int src_rank = (processRank + processCount - 1) % processCount;
    int dest_rank = (processRank + 1) % processCount;
    int current_rank;

    copyVector(replaceVectorXPart, vectorXPart, lineCounts[processRank]);

    for (int i = 0; i < processCount; ++i) {
        current_rank = (processRank + i) % processCount;
        for (int j = 0; j < lineCounts[processRank]; ++j) {
            if (i == 0) {
                matrixAxbPart[j] = -vectorBPart[j];
            }
            for (int k = 0; k < lineCounts[current_rank]; ++k) {
                matrixAxbPart[j] += matrixAPart[j * N + lineOffsets[current_rank] + k] * replaceVectorXPart[k];
            }
        }

        if (i != processCount - 1) {
            MPI_Sendrecv_replace(replaceVectorXPart, lineCounts[0], MPI_DOUBLE, dest_rank, processRank,
                                 src_rank, src_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
}
