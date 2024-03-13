#include <math.h>
#include <mpi/mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 15000
#define EPSILON 1E-7
#define TAU 1E-5
#define MAX_ITERATION_COUNT 100000

void initLineCountsAndOffsets(int* lineCountsPerProcess, int* lineOffsets, int processTotalSize);
void initMatrixAParts(double* matrixAPart, int lineCountInCurrentProcess, int processRank);
void initVectorBParts(double* vectorBPart, int size);
double calculateEuclideanVectorNorm(const double* vector, int size);
void calculate_Axb(const double* matrixAPart, const double* vectorXPart, const double* vectorBPart, double* replaceVectorXPart,
                   double* vectorAxbPart, int* lineCounts, const int* lineOffsets, int processRank, int processTotalSize);
void setVectorXtoFofX(const double* vectorAxb, double* vectorXPart, int partSize);

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
    initLineCountsAndOffsets(lineCounts, lineOffsets, processTotalSize);

    vectorXPart = calloc(lineCounts[processCurrentRank], sizeof(double));
    vectorBPart = malloc(sizeof(double) * lineCounts[processCurrentRank]);
    matrixAPart = malloc(sizeof(double) * lineCounts[processCurrentRank] * N);
    initVectorBParts(vectorBPart, lineCounts[processCurrentRank]);
    initMatrixAParts(matrixAPart, lineCounts[processCurrentRank], lineOffsets[processCurrentRank]);

    vectorBPartNorm = calculateEuclideanVectorNorm(vectorBPart, lineCounts[processCurrentRank]);
    MPI_Reduce(&vectorBPartNorm, &vectorBNorm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (processCurrentRank == 0) {
        printf("All necessary data has been initialized.\nTotal process count: %d\n", processCurrentRank);
        vectorBNorm = sqrt(vectorBNorm);
    }
    vectorAxbPart = malloc(sizeof(double) * lineCounts[processCurrentRank]);
    replaceVectorXParts = malloc(sizeof(double) * lineCounts[0]);


    beginningTime = MPI_Wtime();
    for (iterationsCount = 0; accuracy > EPSILON && iterationsCount < MAX_ITERATION_COUNT; ++iterationsCount) {
        calculate_Axb(matrixAPart, vectorXPart, vectorBPart, replaceVectorXParts, vectorAxbPart, lineCounts, lineOffsets, processCurrentRank,
                      processTotalSize);

        setVectorXtoFofX(vectorAxbPart, vectorXPart, lineCounts[processCurrentRank]);
        vectorAxbNorm = calculateEuclideanVectorNorm(vectorAxbPart, lineCounts[processCurrentRank]);

        //суммирование квадратов норм частей vectorAxb в нулевом процессе для дальнейшего вычисления в нём корня
        MPI_Reduce(&vectorAxbNorm, &accuracy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (processCurrentRank == 0) {
            accuracy = sqrt(accuracy) / vectorBNorm;
        }
        //обновление accuracy на всех процессах для проверки условия
        MPI_Bcast(&accuracy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    endingTime = MPI_Wtime();


    if (processCurrentRank == 0) {
        if (iterationsCount == MAX_ITERATION_COUNT) {
            fprintf(stderr, "Too many iterations\n");
        }else{
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


void initMatrixAParts(double* matrixAPart, int lineCountInCurrentProcess, int processRank) {
    for (int i = 0; i < lineCountInCurrentProcess; ++i) {
        for (int j = 0; j < N; ++j) {
            matrixAPart[i * N + j] = 1;
        }
        matrixAPart[i * N + processRank + i] = 2;
    }
}
void initVectorBParts(double* vectorBPart, int size) {
    for (int i = 0; i < size; ++i) {
        vectorBPart[i] = N + 1;
    }
}
void initLineCountsAndOffsets(int* lineCountsPerProcess, int* lineOffsets, int processTotalSize) {
    int offset = 0;
    for (int i = 0; i < processTotalSize; ++i) {
        lineCountsPerProcess[i] = N / processTotalSize;

        if (i < N % processTotalSize) {
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
void setVectorXtoFofX(const double* vectorAxb, double* vectorXPart, int partSize) {
    for (int i = 0; i < partSize; ++i) {
        vectorXPart[i] -= TAU * vectorAxb[i];
    }
}
void calculate_Axb(const double* matrixAPart, const double* vectorXPart, const double* vectorBPart, double* replaceVectorXPart,
                   double* vectorAxbPart, int* lineCounts, const int* lineOffsets, int processRank, int processTotalSize) {
    int sourceRank = (processRank + processTotalSize - 1) % processTotalSize;
    int destinationRank = (processRank + 1) % processTotalSize;
    int processCurrentRank;

    memcpy(replaceVectorXPart, vectorXPart, lineCounts[processRank]);

    for (int i = 0; i < processTotalSize; ++i) {
        processCurrentRank = (processRank + i) % processTotalSize;
        for (int j = 0; j < lineCounts[processRank]; ++j) {
            if (i == 0) {
                vectorAxbPart[j] = -vectorBPart[j];
            }
            for (int k = 0; k < lineCounts[processCurrentRank]; ++k) {
                vectorAxbPart[j] += matrixAPart[j * N + lineOffsets[processCurrentRank] + k] * replaceVectorXPart[k];
            }
        }

        if (i != processTotalSize - 1) {
            //
            MPI_Sendrecv_replace(replaceVectorXPart, lineCounts[0], MPI_DOUBLE, destinationRank, processRank,
                                 sourceRank, sourceRank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
}
