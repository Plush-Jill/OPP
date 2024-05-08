#include <math.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define DIMS_COUNT 2
#define X 0
#define Y 1

void initDimensions(int dims[2], int processCount);
void initCommunicators(const int dimensions[DIMS_COUNT], MPI_Comm *commGrid, MPI_Comm *commRows, MPI_Comm *commColumns);
void generate_matrix(double *matrix, int column, int leading_row, int leading_column, bool onRows);
void split_A(const double *A, double *A_block, int A_block_size, int n_2, int coords_y, MPI_Comm comm_rows, MPI_Comm comm_columns);
void split_B(const double *B, double *B_block, int B_block_size, int n_2, int aligned_n3, int coords_x, MPI_Comm comm_rows, MPI_Comm comm_columns);
void multiply(const double *A_block, const double *B_block, double *C_block, int A_block_size, int B_block_size, int n_2);
void gather_C(const double *C_block, double *C, int A_block_size, int B_block_size, int aligned_n1, int aligned_n3, int proc_count, MPI_Comm comm_grid);
bool check_C(const double *C, int column, int leading_row, int leading_column, int n_2);

int main(int argc, char **argv) {
    int matrixSize1 = 2000;
    int matrixSize2 = 1500;
    int matrixSize3 = 2500;
    int processRank;
    int processCount;
    int aligned_n1;
    int aligned_n3;
    int matrixABlockSize;
    int matrixBBlockSize;
    int dims[DIMS_COUNT] = {};
    int coords[DIMS_COUNT] = {};
    double beginningTime;
    double endingTime;
    double* A = NULL;
    double* B = NULL;
    double* C = NULL;
    double* A_block = NULL;
    double* B_block = NULL;
    double* C_block = NULL;
    MPI_Comm commGrid;
    MPI_Comm commRows;
    MPI_Comm commColumns;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank);
    MPI_Comm_size(MPI_COMM_WORLD, &processCount);

    initDimensions(dims, processCount);

    initCommunicators(dims, &commGrid, &commRows, &commColumns);

    // Get coordinates of processes
    MPI_Cart_coords(commGrid, processRank, DIMS_COUNT, coords);

    // Set parameters of matrix blocks
    matrixABlockSize = ceil((double)matrixSize1 / dims[X]);
    matrixBBlockSize = ceil((double)matrixSize3 / dims[Y]);
    aligned_n1 = matrixABlockSize * dims[X];
    aligned_n3 = matrixBBlockSize * dims[Y];

    if (coords[X] == 0 && coords[Y] == 0) {
        A = malloc(sizeof(double) * aligned_n1 * matrixSize2);
        B = malloc(sizeof(double) * matrixSize2 * aligned_n3);
        C = malloc(sizeof(double) * aligned_n1 * aligned_n3);

        generate_matrix(A, matrixSize2, matrixSize1, matrixSize2, true);
        generate_matrix(B, aligned_n3, matrixSize2, matrixSize3, false);
    }

    beginningTime = MPI_Wtime();

    A_block = malloc(sizeof(double) * matrixABlockSize * matrixSize2);
    B_block = malloc(sizeof(double) * matrixBBlockSize * matrixSize2);
    C_block = malloc(sizeof(double) * matrixABlockSize * matrixBBlockSize);

    split_A(A, A_block, matrixABlockSize, matrixSize2, coords[Y], commRows, commColumns);
    split_B(B, B_block, matrixBBlockSize, matrixSize2, aligned_n3, coords[X], commRows, commColumns);

    multiply(A_block, B_block, C_block, matrixABlockSize, matrixBBlockSize, matrixSize2);

    gather_C(C_block, C, matrixABlockSize, matrixBBlockSize, aligned_n1, aligned_n3, processCount, commGrid);

    endingTime = MPI_Wtime();

    if (coords[Y] == 0 && coords[X] == 0) {
        printf("Is matrix C correct? - %s\n", check_C(C, aligned_n3, matrixSize1, matrixSize3, matrixSize2) ? "yes" : "no");
        printf("Time: %lf\n", endingTime - beginningTime);

        free(A);
        free(B);
        free(C);
    }

    free(A_block);
    free(B_block);
    free(C_block);
    MPI_Comm_free(&commGrid);
    MPI_Comm_free(&commRows);
    MPI_Comm_free(&commColumns);

    MPI_Finalize();

    return EXIT_SUCCESS;
}

/**
 * @brief Initialize dimensions
 *
 * @param dims
 * @param processCount
 * @param argc
 * @param argv
 */
void initDimensions(int dims[2], int processCount) {
    if (argc < 3) {
        MPI_Dims_create(processCount, DIMS_COUNT, dims);
    } else {
        dims[X] = atoi(argv[1]);
        dims[Y] = atoi(argv[2]);

        if (dims[X] * dims[Y] != processCount) {
            exit(EXIT_FAILURE);
        }
    }
}

/**
 * @brief Initialize communicators
 *
 * @param dims Size of dimensions
 * @param comm_grid Communicator for grids
 * @param comm_rows Comunnicator for rows
 * @param comm_columns Comunnicator for columns
 *
 * @returns comm_grid, comm_rows, comm_columns
 */
void initCommunicators(const int dims[DIMS_COUNT], MPI_Comm *comm_grid, MPI_Comm *comm_rows, MPI_Comm *comm_columns) {
    int reorder = 1;
    int periods[DIMS_COUNT] = {};
    int sub_dims[DIMS_COUNT] = {};

    MPI_Cart_create(MPI_COMM_WORLD, DIMS_COUNT, dims, periods, reorder, comm_grid);

    sub_dims[X] = false;
    sub_dims[Y] = true;
    MPI_Cart_sub(*comm_grid, sub_dims, comm_rows);

    sub_dims[X] = true;
    sub_dims[Y] = false;
    MPI_Cart_sub(*comm_grid, sub_dims, comm_columns);
}

/**
 * @brief Generate matrix that has same numbers on rows or columns
 *
 * @param matrix Pointer to array of row*column size
 * @param columnn Number of column
 * @param leading_row Number of rows in which data will be stored (leading_row <= row)
 * @param leading_column Number of columns in which data will be stored (leading_column <= column)
 * @param onRows true - rows have same numbers, false - columns have same numbers
 *
 * @returns matrix
 */
void generate_matrix(double *matrix, int column, int leading_row, int leading_column, bool onRows) {
    for (int i = 0; i < leading_row; ++i) {
        for (int j = 0; j < leading_column; ++j) {
            matrix[i * column + j] = onRows ? i : j;
        }
    }
}

/**
 * @brief Split matrix A into row blocks between processes
 *
 * @param A Matrix A which is only available in 0 process
 * @param A_block Row block
 * @param A_block_size Number of rows in row block
 * @param n_2 Number of columns of matrix A and number of rows of matrix B
 * @param coords_y Y coordinate of process
 * @param comm_rows Comunnicator for rows
 * @param comm_columns Comunnicator for columns
 *
 * @returns A_block
 */
void split_A(const double *A, double *A_block, int A_block_size, int n_2, int coords_y, MPI_Comm comm_rows, MPI_Comm comm_columns){
    if (coords_y == 0) {
        MPI_Scatter(A, A_block_size * n_2, MPI_DOUBLE, A_block, A_block_size * n_2, MPI_DOUBLE, 0, comm_columns);
    }

    MPI_Bcast(A_block, A_block_size * n_2, MPI_DOUBLE, 0, comm_rows);
}

/**
 * @brief Split matrix B into column blocks between processes
 *
 * @param B Matrix B which is only available in 0 process
 * @param B_block Column block
 * @param B_block_size Number of columns in column block
 * @param n_2 Number of columns of matrix A and number of rows of matrix B
 * @param aligned_n3  Aligned size n_3
 * @param coords_x X coordinate of process
 * @param comm_rows Comunnicator for rows
 * @param comm_columns Comunnicator for columns
 *
 * @returns B_block
 */
void split_B(const double *B, double *B_block, int B_block_size, int n_2, int aligned_n3, int coords_x, MPI_Comm comm_rows, MPI_Comm comm_columns){
    if (coords_x == 0) {
        MPI_Datatype column_not_resized_t;
        MPI_Datatype column_resized_t;

        MPI_Type_vector(n_2, B_block_size, aligned_n3, MPI_DOUBLE, &column_not_resized_t);
        MPI_Type_commit(&column_not_resized_t);

        MPI_Type_create_resized(column_not_resized_t, 0, B_block_size * sizeof(double), &column_resized_t);
        MPI_Type_commit(&column_resized_t);

        MPI_Scatter(B, 1, column_resized_t, B_block, B_block_size * n_2, MPI_DOUBLE, 0, comm_rows);

        MPI_Type_free(&column_not_resized_t);
        MPI_Type_free(&column_resized_t);
    }

    MPI_Bcast(B_block, B_block_size * n_2, MPI_DOUBLE, 0, comm_columns);
}

/**
 * @brief Multiply row block of matrix A and column block of matrix B
 *
 * @param A_block Row block of matrix A
 * @param B_block Column block of matrix B
 * @param C_block Grid block of matrix C
 * @param A_block_size Number of rows in row block
 * @param B_block_size Number of columns in column block
 * @param n_2 Number of columns of matrix A and number of rows of matrix B
 *
 * @returns C_block
 */
void multiply(const double *A_block, const double *B_block, double *C_block, int A_block_size, int B_block_size, int n_2){
    for (int i = 0; i < A_block_size; ++i) {
        for (int j = 0; j < B_block_size; ++j) {
            C_block[i * B_block_size + j] = 0;
        }
    }

    for (int i = 0; i < A_block_size; ++i) {
        for (int j = 0; j < n_2; ++j) {
            for (int k = 0; k < B_block_size; ++k) {
                C_block[i * B_block_size + k] += A_block[i * n_2 + j] * B_block[j * B_block_size + k];
            }
        }
    }
}

/**
 * @brief Gather matrix C of grid blocks
 *
 * @param C_block Grid block of matrix C
 * @param C Matrix C
 * @param A_block_size Number of rows in row block
 * @param B_block_size Number of columns in column block
 * @param aligned_n1 Aligned number of number of rows of matrix A
 * @param aligned_n3 Aligned number of columns of matrix B
 * @param proc_rank Rank of current process
 * @param comm_grid Communicator for grids
 *
 * @returns C
 */
void gather_C(const double *C_block, double *C, int A_block_size, int B_block_size, int aligned_n1, int aligned_n3, int proc_count, MPI_Comm comm_grid){
    MPI_Datatype not_resized_recv_t;
    MPI_Datatype resized_recv_t;

    int dims_x = aligned_n1 / A_block_size;
    int dims_y = aligned_n3 / B_block_size;
    int *recv_counts = malloc(sizeof(int) * proc_count);
    int *displs = malloc(sizeof(int) * proc_count);

    MPI_Type_vector(A_block_size, B_block_size, aligned_n3, MPI_DOUBLE, &not_resized_recv_t);
    MPI_Type_commit(&not_resized_recv_t);

    MPI_Type_create_resized(not_resized_recv_t, 0, B_block_size * sizeof(double), &resized_recv_t);
    MPI_Type_commit(&resized_recv_t);

    for (int i = 0; i < dims_x; ++i) {
        for (int j = 0; j < dims_y; ++j) {
            recv_counts[i * dims_y + j] = 1;
            displs[i * dims_y + j] = j + i * dims_y * A_block_size;
        }
    }

    MPI_Gatherv(C_block, A_block_size * B_block_size, MPI_DOUBLE, C, recv_counts, displs, resized_recv_t, 0, comm_grid);

    MPI_Type_free(&not_resized_recv_t);
    MPI_Type_free(&resized_recv_t);
    free(recv_counts);
    free(displs);
}

/**
 * @brief Сheck result of multiplying matrices A and B stored in matrix C for correctness
 *
 * @param C matrix C storing rsult of multiplying matrices A and B
 * @param column Number of column
 * @param leading_row Number of rows in which data will be stored (leading_row <= row)
 * @param leading_column Number of columns in which data will be stored (leading_column <= column)
 * @param n_2 Number of columns of matrix A and number of rows of matrix B
 *
 * @return true - Result of multiplying matrices A and B stored in matrix C is correct,
 * @return false - Result of multiplying matrices A and B stored in matrix C is incorrect
 */
bool check_C(const double *C, int column, int leading_row, int leading_column, int n_2) {
    for (int i = 0; i < leading_row; ++i) {
        for (int j = 0; j < leading_column; ++j) {
            if (C[i * column + j] != (double) i * j * n_2) {
                printf("(%d, %d)\n", i, j);
                printf("%lf != %lf\n", C[i * column + j], (double) i * j * n_2);
                return false;
            }
        }
    }

    return true;
}