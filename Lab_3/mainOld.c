#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <mpi/mpi.h>

#define DIMS_COUNT 2
#define X 0
#define Y 1

#define N1 4320
#define N2 9000
#define N3 2280

static int rank, sizeProccess;

void breakProgramm() {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    exit(-1);
}

void printMatrix(double *matrix, int n1, int n2) {
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            printf("%f ", matrix[i * n2 + j]);
        }
        printf("\n");
    }
    printf("\n");
}

/*----------*/

void fillMatrix(double *matrix, int n1, int n2) {
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            matrix[i * n2 + j] = i + j;
        }
    }
}

void initCommunicators(const int dims[DIMS_COUNT], MPI_Comm *commGrid, MPI_Comm *commRows, MPI_Comm *commColumns) {
    int reorder = 1;
    int periods[DIMS_COUNT] = {};
    int subDims[DIMS_COUNT] = {};

    MPI_Cart_create(MPI_COMM_WORLD, DIMS_COUNT, dims, periods, reorder, commGrid);

    subDims[X] = false;
    subDims[Y] = true;
    MPI_Cart_sub(*commGrid, subDims, commRows);

    subDims[X] = true;
    subDims[Y] = false;
    MPI_Cart_sub(*commGrid, subDims, commColumns);
}

void splitA(double *matrix1, double *matrix1Block, int matrix1BlockSize,
            int n2, int coordsY, MPI_Comm commRows, MPI_Comm commColumns) {

    if (coordsY == 0) {

        MPI_Scatter(matrix1, matrix1BlockSize * n2, MPI_DOUBLE, matrix1Block, matrix1BlockSize * n2, MPI_DOUBLE, 0, commColumns);

    }

    MPI_Bcast(matrix1Block, matrix1BlockSize * n2, MPI_DOUBLE, 0, commRows);
}

void splitB(double *matrix2, double *matrix2Block, int matrix2BlockSize,
            int n2, int alignedN3, int coordsX, MPI_Comm commRows, MPI_Comm commColumns) {

    if (coordsX == 0) {

        MPI_Datatype columnNotResized;
        MPI_Datatype columnResized;

        MPI_Type_vector(N2, matrix2BlockSize, alignedN3, MPI_DOUBLE, &columnNotResized);
        MPI_Type_commit(&columnNotResized);

        MPI_Type_create_resized(columnNotResized, 0, matrix2BlockSize * sizeof(double), &columnResized);
        MPI_Type_commit(&columnResized);

        MPI_Scatter(matrix2, 1, columnResized, matrix2Block, matrix2BlockSize * N2, MPI_DOUBLE, 0, commRows);

        MPI_Type_free(&columnNotResized);
        MPI_Type_free(&columnResized);

    }

    MPI_Bcast(matrix2Block, matrix2BlockSize * n2, MPI_DOUBLE, 0, commColumns);
}

void multiply(double *matrix1Block, double *matrix2Block, double *matrix3Block,
              int matrix1BlockSize, int matrix2BlockSize, int n2) {

    memset(matrix3Block, 0, matrix1BlockSize * matrix2BlockSize * sizeof(double));

    for (int i = 0; i < matrix1BlockSize; ++i) {
        for (int j = 0; j < n2; ++j) {
            for (int k = 0; k < matrix2BlockSize; ++k) {
                matrix3Block[i * matrix2BlockSize + k] += matrix1Block[i * n2 + j] * matrix2Block[j * matrix2BlockSize + k];
            }
        }
    }
}

void gatherV(double *matrix3Block, double *matrix3, int matrix1BlockSize, int matrix2BlockSize,
             int alignedN1, int alignedN3, MPI_Comm commGrid) {

    MPI_Datatype notResizedRecv;
    MPI_Datatype resizedRecv;

    int dimsX = alignedN1 / matrix1BlockSize;
    int dimsY = alignedN3 / matrix2BlockSize;
    int *recvCounts = malloc(sizeProccess * sizeof(int));
    int *displs = malloc(sizeProccess * sizeof(int));

    MPI_Type_vector(matrix1BlockSize, matrix2BlockSize, alignedN3, MPI_DOUBLE, &notResizedRecv);
    MPI_Type_commit(&notResizedRecv);

    MPI_Type_create_resized(notResizedRecv, 0, matrix2BlockSize * sizeof(double), &resizedRecv);
    MPI_Type_commit(&resizedRecv);

    for (int i = 0; i < dimsX; ++i) {
        for (int j = 0; j < dimsY; ++j) {
            recvCounts[i * dimsY + j] = 1;
            displs[i * dimsY + j] = j + i * dimsY * matrix1BlockSize;
        }
    }
    MPI_Gatherv(matrix3Block, matrix1BlockSize * matrix2BlockSize, MPI_DOUBLE, matrix3, recvCounts, displs, resizedRecv, 0, commGrid);

    MPI_Type_free(&notResizedRecv);
    MPI_Type_free(&resizedRecv);
    free(recvCounts);
    free(displs);
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &sizeProccess);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int dimsX = (int)sqrt(sizeProccess);
    int dimsY = (int)sqrt(sizeProccess);

    int dims[DIMS_COUNT] = {dimsX, dimsY};
    MPI_Dims_create(sizeProccess, DIMS_COUNT, dims);

    MPI_Comm commGrid;
    MPI_Comm commRows;
    MPI_Comm commColumns;
    initCommunicators(dims, &commGrid, &commRows, &commColumns);

    int coords[DIMS_COUNT] = {};
    MPI_Cart_coords(commGrid, rank, DIMS_COUNT, coords);

    double *matrix1 = NULL;
    double *matrix2 = NULL;
    double *matrix3 = NULL;

    int matrix1BlockSize = N1 / dims[X];
    int matrix2BlockSize = N3 / dims[Y];
    int alignedN1 = matrix1BlockSize * dims[X];
    int alignedN3 = matrix2BlockSize * dims[Y];
    if (coords[X] == 0 && coords[Y] == 0) {

        matrix1 = malloc(alignedN1 * N2 * sizeof(double));
        matrix2 = malloc(N2 * alignedN3 * sizeof(double));
        matrix3 = malloc(alignedN1 * alignedN3 * sizeof(double));

        fillMatrix(matrix1, alignedN1, N2);
        fillMatrix(matrix2, N2, alignedN3);

    }

    // if (rank == 0) {

    //     printMatrix(matrix1, alignedN1, N2);
    //     printMatrix(matrix2, N2, alignedN3);

    // }

    double startTime = MPI_Wtime();

    double *matrix1Block = malloc(matrix1BlockSize * N2 * sizeof(double));
    double *matrix2Block = malloc(matrix2BlockSize * N2 * sizeof(double));
    double *matrix3Block = malloc(matrix1BlockSize * matrix2BlockSize * sizeof(double));

    splitA(matrix1, matrix1Block, matrix1BlockSize, N2, coords[Y], commRows, commColumns);
    splitB(matrix2, matrix2Block, matrix2BlockSize, N2, alignedN3, coords[X], commRows, commColumns);

    multiply(matrix1Block, matrix2Block, matrix3Block, matrix1BlockSize, matrix2BlockSize, N2);

    gatherV(matrix3Block, matrix3, matrix1BlockSize, matrix2BlockSize, alignedN1, alignedN3, commGrid);

    double finishTime = MPI_Wtime();

    // if (rank == 0) {
    //     printMatrix(matrix3, alignedN1, alignedN3);
    // }

    if (rank == 0) {

        printf("%f\n", finishTime - startTime);

        free(matrix1);
        free(matrix2);
        free(matrix3);

    }

    free(matrix1Block);
    free(matrix2Block);
    free(matrix3Block);
    MPI_Comm_free(&commGrid);
    MPI_Comm_free(&commRows);
    MPI_Comm_free(&commColumns);

    MPI_Finalize();
    return 0;
}