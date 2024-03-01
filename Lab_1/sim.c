#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi/mpi.h>

#define PI 3.14159265358979323846
#define N 500
#define epsilon (double)0.00001
#define tao (double)0.0003

static int rank, sizeProcess;
void printMatrix(double *matrix) {
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            printf("%f ", matrix[i * N + j]);
        }
        printf("\n");
    }
}
void printVector(double *vector) {
    for (size_t i = 0; i < N; ++i) {
        printf("%f ", vector[i]);
    }
    printf("\n");
}
void setZeroVector(double *vector) {
    memset(vector, 0, N * sizeof(double));
}
void subVector(double *vector1, double *vector2) {
    for (size_t i = 0; i < N; ++i) {
        vector1[i] -= vector2[i];
    }
}
void multiplyMatrixWithVector(double *pieceMatrix, double *inputVector, double *outputVector,
                              double *vectorBuffer, MPI_Status st) {

    int vectorQuantity = N / sizeProcess;
    if ((N % sizeProcess != 0) && (N % sizeProcess >= rank + 1)) {
        vectorQuantity++;
    }

    if (rank != 0) {
        inputVector = vectorBuffer;
    }
    MPI_Bcast(inputVector, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        for (size_t i = 0; i < vectorQuantity; ++i) {
            for (size_t j = 0; j < N; ++j) {
                outputVector[i] += pieceMatrix[i * N + j] * inputVector[j];
            }
        }
    }
    if (rank != 0) {
        double res = 0;
        for (size_t i = 0; i < vectorQuantity; ++i) {
            res = 0;
            for (size_t j = 0; j < N; ++j) {
                res += pieceMatrix[i * N + j] * inputVector[j];
            }
            MPI_Send(&res, 1, MPI_DOUBLE, 0, 1992, MPI_COMM_WORLD);
        }
    }

    if (rank == 0) {
        double res = 0;
        size_t posInOutputVector = vectorQuantity;
        for (size_t i = 1; i < sizeProcess; ++i) {
            int vectorQuantityInAnotherProcess = N / sizeProcess;
            if ((N % sizeProcess != 0) && (N % sizeProcess >= i + 1)) {
                vectorQuantityInAnotherProcess++;
            }

            for (size_t j = 0; j < vectorQuantityInAnotherProcess; ++j) {
                MPI_Recv(&res, 1, MPI_DOUBLE, i, 1992, MPI_COMM_WORLD, &st);
                outputVector[posInOutputVector] = res;
                ++posInOutputVector;
            }
        }
    }
}
void F_ofXn(double* vectorX, const double* vectorAmult_XpowN__minusB){
    if (rank == 0) {
        for (size_t i = 0; i < N; ++i) {
            vectorX[i] = vectorX[i] - (tao * vectorAmult_XpowN__minusB[i]);
        }
    }
}
double G_ofXn(double* pieceMatrix, double* vectorX, double* vectorB, double* vectorAmult_XpowN__minusB, double* vectorBuffer, MPI_Status* st) {
    subVector(vectorAmult_XpowN__minusB, vectorB);
    double numerator = 0, denominator = 0;
    for (size_t i = 0; i < N; ++i) {
        double x = vectorAmult_XpowN__minusB[i];
        numerator += (x * x);
    }
    numerator = sqrt(numerator);
    for (size_t i = 0; i < N; ++i) {
        double x = vectorB[i];
        denominator += (x * x);
    }
    denominator = sqrt(denominator);

    return (numerator / denominator);
}



int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Status st;
    MPI_Comm_size(MPI_COMM_WORLD, &sizeProcess);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double *pieceMatrix = NULL;
    double *vectorU = NULL;
    double *vectorX = NULL;
    double *vectorB = NULL;
    double *vectorBuffer = NULL;
    double *vectorAmult_XpowN__minusB = NULL;


    int vectorQuantity = N / sizeProcess;
    if ((N % sizeProcess != 0) && (N % sizeProcess >= rank + 1)) {
        vectorQuantity++;
    }
    int vectorQuantityInPrevProcess = 0;
    for (size_t i = 0; i < rank; ++i) {
        vectorQuantityInPrevProcess += N / sizeProcess;
        if ((N % sizeProcess != 0) && (N % sizeProcess >= i + 1)) {
            vectorQuantityInPrevProcess++;
        }
    }
    pieceMatrix = calloc(vectorQuantity * N, sizeof(double));
    for (size_t i = 0; i < vectorQuantity; ++i) {
        for (size_t j = 0; j < N; ++j) {
            pieceMatrix[i * N + j] = 1;
            if (j == i + vectorQuantityInPrevProcess) {
                pieceMatrix[i * N + j] = 2;
            }
        }
    }


    if (rank == 0) {
        printf("MPI\n");
        vectorU = calloc(N, sizeof(double));
        for (size_t i = 0; i < N; ++i) {
            vectorU[i] = sin(2 * PI * (i + 1) / N);
        }

        vectorX = calloc(N, sizeof(double));
        vectorB = calloc(N, sizeof(double));
        vectorAmult_XpowN__minusB = calloc(N, sizeof(double));
        setZeroVector(vectorB);
    }else{
        vectorBuffer = calloc(N, sizeof(double));
    }


    multiplyMatrixWithVector(pieceMatrix, vectorU, vectorB, vectorBuffer, st);

    double startTime = 0;
    if (rank == 0) {
        startTime = MPI_Wtime();
    }

    int isComplete = 0;
    for(size_t k = 0; 1; ++k) {
        if (rank == 0) {
            setZeroVector(vectorAmult_XpowN__minusB);
        }

        multiplyMatrixWithVector(pieceMatrix, vectorX, vectorAmult_XpowN__minusB, vectorBuffer, st);

        if (rank == 0) {
            if (G_ofXn(pieceMatrix, vectorX, vectorB, vectorAmult_XpowN__minusB, vectorBuffer, &st) < epsilon){
                isComplete = 1;
            }
            for (size_t i = 1; i < sizeProcess; ++i) {
                MPI_Send(&isComplete, 1, MPI_INT, i, 199, MPI_COMM_WORLD);
            }
        }else{
            MPI_Recv(&isComplete, 1, MPI_INT, 0, 199, MPI_COMM_WORLD, &st);
        }

        if (isComplete) {
            if (rank == 0) {
                printf("%ld\n", k);
            }
            break;
        }
        F_ofXn(vectorX, vectorAmult_XpowN__minusB);
    }

    if (rank == 0) {
        double endTime = MPI_Wtime();
        printf("%f\n", endTime - startTime);
    }
    free(pieceMatrix);
    free(vectorBuffer);
    free(vectorU);
    free(vectorX);
    free(vectorB);
    free(vectorAmult_XpowN__minusB);

    MPI_Finalize();
    return 0;
}
