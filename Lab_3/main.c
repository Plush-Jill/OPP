#include <mpi/mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define DIMENSIONS_COUNT 2
#define X 0
#define Y 1


void initCommunicators(const int dimensions[DIMENSIONS_COUNT], MPI_Comm *commGrid, MPI_Comm *commRows, MPI_Comm *commColumns);
void initMatrix(double *matrix, int columns, int rows);
void splitMatrixA(const double *matrixA, double *matrixABlock, int matrixABlockSize, int n2, int coordsY, MPI_Comm commRows, MPI_Comm commColumns);
void splitMatrixB(const double *matrixB, double *matrixBBlock, int matrixBBlockSize, int n2, int n3, int coordsX, MPI_Comm commRows, MPI_Comm commColumns);
void multiply(const double *matrixABlock, const double *matrixBBlock, double *matrixCBlock, int matrixABlockSize,
              int matrixBBlockSize, int n2, int processRank);
void gatherMatrixC(const double *matrixCBlock, double *matrixC, int matrixABlockSize, int matrixBBlockSize, int n1, int n3,
              int processCount, MPI_Comm commGrid, int processRank);
bool checkMatrixC(const double *matrixC, int columns, int rows, int n2);
void printMatrix(const double* matrix, int rows, int columns);



int main(int argc, char **argv) {
    int n1 = 144 * 20;
    int n2 = 144 * 24;
    int n3 = 144 * 25;
    int processRank;
    int processCount;
    int matrixABlockSize;
    int matrixBBlockSize;
    int dimensions[DIMENSIONS_COUNT] = {};
    int coords[DIMENSIONS_COUNT] = {};
    double beginningTime;
    double endingTime;
    double* matrixA = NULL;
    double* matrixB = NULL;
    double* matrixC = NULL;
    double* matrixABlock = NULL;
    double* matrixBBlock = NULL;
    double* matrixCBlock = NULL;
    MPI_Comm commGrid;
    MPI_Comm commRows;
    MPI_Comm commColumns;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank);
    MPI_Comm_size(MPI_COMM_WORLD, &processCount);
    if (processRank == 0) {
        fprintf(stdout, "Process count: %d\n", processCount);
    }

    MPI_Dims_create(processCount, DIMENSIONS_COUNT, dimensions);

    initCommunicators(dimensions, &commGrid, &commRows, &commColumns);
    // Получение координат текущего процесса в comm.
    MPI_Cart_coords(commGrid, processRank, DIMENSIONS_COUNT, coords);

    matrixABlockSize = n1 / dimensions[X];
    matrixBBlockSize = n3 / dimensions[Y];

    if (coords[X] == 0 && coords[Y] == 0) {
        if (processRank != 0) {
            fprintf(stderr, "(0, 0) isn't 0'th process.");
        }

        matrixA = malloc(sizeof(double) * n1 * n2);
        matrixB = malloc(sizeof(double) * n2 * n3);
        matrixC = malloc(sizeof(double) * n1 * n3);

        initMatrix(matrixA, n1, n2);
        initMatrix(matrixB, n2, n3);
    }
    /*if (processRank == 0) {
        fprintf(stdout, "Coordinates:\n");
    }*/
    //MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < processCount; ++i) {
        if (processRank == i) {
            fprintf(stdout, "Process %d: (%d, %d).\n", processRank, coords[X], coords[Y]);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    matrixABlock = malloc(sizeof(double) * matrixABlockSize * n2);
    matrixBBlock = malloc(sizeof(double) * matrixBBlockSize * n2);
    matrixCBlock = calloc(matrixABlockSize * matrixBBlockSize, sizeof(double));

    beginningTime = MPI_Wtime();

    splitMatrixA(matrixA, matrixABlock, matrixABlockSize, n2, coords[Y], commRows, commColumns);
    splitMatrixB(matrixB, matrixBBlock, matrixBBlockSize, n2, n3, coords[X], commRows, commColumns);

    multiply(matrixABlock, matrixBBlock, matrixCBlock, matrixABlockSize, matrixBBlockSize, n2, processRank);

    gatherMatrixC(matrixCBlock, matrixC, matrixABlockSize, matrixBBlockSize, n1, n3, processCount, commGrid, processRank);

    endingTime = MPI_Wtime();

    if (coords[Y] == 0 && coords[X] == 0) {
        printf("MatrixC is %s\n", checkMatrixC(matrixC, n3, n1, n2) ? "correct" : "incorrect");
        printf("Time: %lf\n", endingTime - beginningTime);

        free(matrixA);
        free(matrixB);
        free(matrixC);
    }

    free(matrixABlock);
    free(matrixBBlock);
    free(matrixCBlock);
    MPI_Comm_free(&commGrid);
    MPI_Comm_free(&commRows);
    MPI_Comm_free(&commColumns);

    MPI_Finalize();

    return EXIT_SUCCESS;
}


void initCommunicators(const int dimensions[DIMENSIONS_COUNT], MPI_Comm* commGrid, MPI_Comm* commRows, MPI_Comm* commColumns) {
    int periods[DIMENSIONS_COUNT] = {};
    int subDimensions[DIMENSIONS_COUNT] = {};

    MPI_Cart_create(MPI_COMM_WORLD, DIMENSIONS_COUNT, dimensions, periods, 1, commGrid);

    subDimensions[X] = false;
    subDimensions[Y] = true;
    MPI_Cart_sub(*commGrid, subDimensions, commRows);

    subDimensions[X] = true;
    subDimensions[Y] = false;
    MPI_Cart_sub(*commGrid, subDimensions, commColumns);
}
void initMatrix(double *matrix, int rows, int columns) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            matrix[i * columns + j] = i;
        }
    }
}
void splitMatrixA(const double* matrixA, double* matrixABlock, int matrixABlockSize, int n2, int coordsY, MPI_Comm commRows, MPI_Comm commColumns){
    if (coordsY == 0) {
        MPI_Scatter(matrixA, matrixABlockSize * n2, MPI_DOUBLE, matrixABlock, matrixABlockSize * n2, MPI_DOUBLE, 0, commColumns);
    }

    MPI_Bcast(matrixABlock, matrixABlockSize * n2, MPI_DOUBLE, 0, commRows);
}
void splitMatrixB(const double* matrixB, double* matrixBBlock, int matrixBBlockSize, int n2, int n3, int coordsX, MPI_Comm commRows, MPI_Comm commColumns){
    if (coordsX == 0) {
        MPI_Datatype columnTypeNotResized;
        MPI_Datatype columnTypeResized;

        MPI_Type_vector(n2, matrixBBlockSize, n3, MPI_DOUBLE, &columnTypeNotResized);
        MPI_Type_commit(&columnTypeNotResized);

        MPI_Type_create_resized(columnTypeNotResized, 0, matrixBBlockSize * sizeof(double), &columnTypeResized);
        MPI_Type_commit(&columnTypeResized);

        MPI_Scatter(matrixB, 1, columnTypeResized, matrixBBlock, matrixBBlockSize * n2, MPI_DOUBLE, 0, commRows);

        MPI_Type_free(&columnTypeNotResized);
        MPI_Type_free(&columnTypeResized);
    }

    MPI_Bcast(matrixBBlock, matrixBBlockSize * n2, MPI_DOUBLE, 0, commColumns);
}
void multiply(const double *matrixABlock, const double *matrixBBlock, double *matrixCBlock, int matrixABlockSize,
              int matrixBBlockSize, int n2, int processRank) {
    long long sum = 0;
    for (int i = 0; i < matrixABlockSize; ++i) {
        for (int j = 0; j < n2; ++j) {
            for (int k = 0; k < matrixBBlockSize; ++k) {
                matrixCBlock[i * matrixBBlockSize + k] += matrixABlock[i * n2 + j] * matrixBBlock[j * matrixBBlockSize + k];
                ++sum;
            }
        }
    }
    /*for (int i = 0; i < 16; ++i){
        if (i == processRank) {
            fprintf(stdout, "process %d sum: %lld\n", processRank, sum);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }*/
}
void gatherMatrixC(const double *matrixCBlock, double *matrixC, int matrixABlockSize, int matrixBBlockSize, int n1, int n3,
              int processCount, MPI_Comm commGrid, int processRank) {
    MPI_Datatype receiveTypeNotResized;
    MPI_Datatype receiveTypeResized;

    int horizontalLinesCount = n1 / matrixABlockSize;
    int verticalLinesCount = n3 / matrixBBlockSize;
    int* receiveCounts = malloc(sizeof(int) * processCount);
    int* offsets = malloc(sizeof(int) * processCount);

    MPI_Type_vector(matrixABlockSize, matrixBBlockSize, n3, MPI_DOUBLE, &receiveTypeNotResized);
    MPI_Type_commit(&receiveTypeNotResized);

    MPI_Type_create_resized(receiveTypeNotResized, 0, matrixBBlockSize * sizeof(double), &receiveTypeResized);
    MPI_Type_commit(&receiveTypeResized);

    for (int i = 0; i < horizontalLinesCount; ++i) {
        for (int j = 0; j < verticalLinesCount; ++j) {
            receiveCounts[i * verticalLinesCount + j] = 1;
            offsets[i * verticalLinesCount + j] = i * verticalLinesCount * matrixABlockSize + j;
        }
    }
    for (int i = 0; i < processCount; ++i){
        if (i == processRank) {
            fprintf(stdout, "process %d: ", processRank);
            for (int j = 0; j < processCount; ++j) {
                fprintf(stdout, "%d ", offsets[j]);
            }
            fprintf(stdout, "\n");
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
    //MPI_Gather(matrixCBlock, matrixABlockSize * matrixBBlockSize, MPI_DOUBLE, matrixC, 1, receiveTypeResized, 0, commGrid);
    MPI_Gatherv(matrixCBlock, matrixABlockSize * matrixBBlockSize, MPI_DOUBLE, matrixC, receiveCounts, offsets, receiveTypeResized, 0, commGrid);

    MPI_Type_free(&receiveTypeNotResized);
    MPI_Type_free(&receiveTypeResized);
    free(receiveCounts);
    free(offsets);
}
bool checkMatrixC(const double *matrixC, int columns, int rows, int n2) {
    long long n2Sum = 0;
    for (int i = 0; i < n2; ++i){
        n2Sum += i;
    }
    long long currentSum = 0;
    for (int i = 0; i < rows; ++i){
        for (int j = 0; j < columns; ++j) {
            if (matrixC[i * columns + j] != currentSum) {
                return false;
            }
        }
        currentSum += n2Sum;
    }

    return true;
}
void printMatrix(const double* matrix, int rows, int columns){
    fprintf(stdout, "printing %dx%d matrix..\n", rows, columns);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            fprintf(stdout, "%.0f ",matrix[i * columns + j]);
        }
        fprintf(stdout, "|\n");
    }
}
