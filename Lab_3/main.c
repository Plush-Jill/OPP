#include <math.h>
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
    // Get coordinates of processes
    MPI_Cart_coords(commGrid, processRank, DIMENSIONS_COUNT, coords);

    // Set parameters of matrix blocks
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
    beginningTime = MPI_Wtime();

    matrixABlock = malloc(sizeof(double) * matrixABlockSize * n2);
    matrixBBlock = malloc(sizeof(double) * matrixBBlockSize * n2);
    matrixCBlock = calloc(matrixABlockSize * matrixBBlockSize, sizeof(double));

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

/**
 * @brief Initialize communicators
 *
 * @param dimensions Size of dimensions
 * @param commGrid Communicator for grids
 * @param commRows Communicator for rows
 * @param commColumns Communicator for columns
 *
 * @returns comm_grid, comm_rows, comm_columns
 */
void initCommunicators(const int dimensions[DIMENSIONS_COUNT], MPI_Comm* commGrid, MPI_Comm* commRows, MPI_Comm* commColumns) {
    int reorder = 1;
    int periods[DIMENSIONS_COUNT] = {};
    int subDimensions[DIMENSIONS_COUNT] = {};

    MPI_Cart_create(MPI_COMM_WORLD, DIMENSIONS_COUNT, dimensions, periods, reorder, commGrid);

    subDimensions[X] = false;
    subDimensions[Y] = true;
    MPI_Cart_sub(*commGrid, subDimensions, commRows);

    subDimensions[X] = true;
    subDimensions[Y] = false;
    MPI_Cart_sub(*commGrid, subDimensions, commColumns);
}

/**
 * @brief Generate matrix that has same numbers on rows or columns
 *
 * @param matrix Pointer to array of row*column size
 * @param columns Number of column
 * @param rows Number of rows in which data will be stored (leadingRow <= row)
 * @param leadingColumn Number of columns in which data will be stored (leadingColumn <= column)
 * @param onRows true - rows have same numbers, false - columns have same numbers
 *
 * @returns matrix
 */
void initMatrix(double *matrix, int rows, int columns) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            matrix[i * columns + j] = i;
        }
    }
}

/**
 * @brief Split matrix A into row blocks between processes
 *
 * @param matrixA Matrix A which is only available in 0 process
 * @param matrixABlock Row block
 * @param matrixABlockSize Number of rows in row block
 * @param n2 Number of columns of matrix A and number of rows of matrix B
 * @param coordsY Y coordinate of process
 * @param commRows Comunnicator for rows
 * @param commColumns Comunnicator for columns
 *
 * @returns A_block
 */
void splitMatrixA(const double* matrixA, double* matrixABlock, int matrixABlockSize, int n2, int coordsY, MPI_Comm commRows, MPI_Comm commColumns){
    if (coordsY == 0) {
        MPI_Scatter(matrixA, matrixABlockSize * n2, MPI_DOUBLE, matrixABlock, matrixABlockSize * n2, MPI_DOUBLE, 0, commColumns);
    }

    MPI_Bcast(matrixABlock, matrixABlockSize * n2, MPI_DOUBLE, 0, commRows);
}

/**
 * @brief Split matrix B into column blocks between processes
 *
 * @param matrixB Matrix B which is only available in 0 process
 * @param matrixBBlock Column block
 * @param matrixBBlockSize Number of columns in column block
 * @param n2 Number of columns of matrix A and number of rows of matrix B
 * @param n3  Aligned size n_3
 * @param coordsX X coordinate of process
 * @param commRows Comunnicator for rows
 * @param commColumns Comunnicator for columns
 *
 * @returns B_block
 */
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

/**
 * @brief Multiply row block of matrix A and column block of matrix B
 *
 * @param matrixABlock Row block of matrix A
 * @param matrixBBlock Column block of matrix B
 * @param matrixCBlock Grid block of matrix C
 * @param matrixABlockSize Number of rows in row block
 * @param matrixBBlockSize Number of columns in column block
 * @param n2 Number of columns of matrix A and number of rows of matrix B
 *
 * @returns C_block
 */
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

/**
 * @brief Gather matrix C of grid blocks
 *
 * @param matrixCBlock Grid block of matrix C
 * @param matrixC Matrix C
 * @param matrixABlockSize Number of rows in row block
 * @param matrixBBlockSize Number of columns in column block
 * @param n1 Aligned number of number of rows of matrix A
 * @param n3 Aligned number of columns of matrix B
 * @param proc_rank Rank of current process
 * @param commGrid Communicator for grids
 *
 * @returns C
 */
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

/**
 * @brief Сheck result of multiplying matrices A and B stored in matrix C for correctness
 *
 * @param matrixC matrix C storing rsult of multiplying matrices A and B
 * @param columns Number of column
 * @param rows Number of rows in which data will be stored (leading_row <= row)
 * @param leading_column Number of columns in which data will be stored (leading_column <= column)
 * @param n_2 Number of columns of matrix A and number of rows of matrix B
 *
 * @return true - Result of multiplying matrices A and B stored in matrix C is correct,
 * @return false - Result of multiplying matrices A and B stored in matrix C is incorrect
 */
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
