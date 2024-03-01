#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi/mpi.h>

#define PI 3.14159265358979323846
#define N 500

static int rank, sizeProccess;
const double epsilon = 0.00001;
const double tao = 0.0003;


void printMatrix(double *matrix) {
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            printf("%f ", matrix[i * N + j]);
        }
        printf("\n");
    }
}
void breakProgramm() {
    MPI_Barrier(MPI_COMM_WORLD);
    exit(-1);
}
void printVector(double *vector, size_t vectorSize) {
    for (size_t i = 0; i < vectorSize; ++i) {
        printf("%f ", vector[i]);
    }
    printf("\n");
}
void printVectorv2(double *vector, size_t vectorSize) {
    for (size_t i = 0; i < sizeProccess; ++i) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (i == rank) {
            // printf("rank %d: ", rank);
            for (size_t j = 0; j < vectorSize; ++j) {
                printf("%f ", vector[j]);
            }
            // printf("\n");
        }
    }
    if (rank == sizeProccess - 1) {
        printf("\n");
        printf("\n");
    }
}
void setZeroVector(double *vector, size_t vectorSize) {
    memset(vector, 0, vectorSize * sizeof(double));
}
void subVector(double *vector1, double *vector2, size_t sizeVector) {
    for (size_t i = 0; i < sizeVector; ++i) {
        vector1[i] -= vector2[i];
    }
}
double getNorm(double *vector, size_t sizeVector) {
    double sum = 0;
    for (size_t i = 0; i < sizeVector; ++i) {
        double a = vector[i];
        sum += (a * a);
    }

    double res = 0;
    MPI_Allreduce(&sum, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sqrt(res);
}
void mulMatrixVector(double *pieceVector, double *circleVector, double *outputVector,
                     size_t vectorSizeInCurrentProcess, size_t shiftSize, size_t sumSizeVectorInPrevProcesses) {

    MPI_Request req[2];
    MPI_Status st;

    size_t indexInVector = sumSizeVectorInPrevProcesses; //позиция в строке в которой работаем
    size_t numberBlockInVector = rank; // номер блока в котором работаем
    for (size_t i = 0; i < sizeProccess; ++i) {
        indexInVector %= N;

        for (size_t j = 0; j < vectorSizeInCurrentProcess; ++j) {
            for (size_t k = 0; k < shiftSize; ++k) {
                outputVector[j] += pieceVector[(j * N + k + indexInVector) % (N * vectorSizeInCurrentProcess)] * circleVector[k];
            }
        }

        size_t blockSize = N / sizeProccess;
        if ((N % sizeProccess != 0) && (N % sizeProccess >= numberBlockInVector + 1)) {
            blockSize++;
        }

        if ((N % sizeProccess != 0) && (shiftSize > blockSize)) {
            indexInVector--;
        }
        indexInVector += shiftSize;
        numberBlockInVector++;
        numberBlockInVector %= sizeProccess;

        /*сдвиг циклический*/
        if (sizeProccess > 1) {
            MPI_Isend(circleVector, shiftSize, MPI_DOUBLE, (rank - 1 + sizeProccess) % sizeProccess, 12345, MPI_COMM_WORLD, &req[0]);
            MPI_Irecv(circleVector, shiftSize, MPI_DOUBLE, (rank + 1	  	 	   ) % sizeProccess, 12345, MPI_COMM_WORLD, &req[1]);
            MPI_Waitall(2, req, &st);
        }
    }
}
void F_ofXn(double* vectorX, double* vectorAxn_b, size_t vectorSizeInCurrentProcess){
    for (size_t i = 0; i < vectorSizeInCurrentProcess; ++i){
        vectorX[i] = vectorX[i] - (tao * vectorAxn_b[i]);
    }
}


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Status st;
    MPI_Comm_size(MPI_COMM_WORLD, &sizeProccess);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double *pieceVectorOfMatrix = NULL;
    size_t vectorSizeInCurrentProcess = N / sizeProccess;
    if ((N % sizeProccess != 0) && (N % sizeProccess >= rank + 1)) {
        vectorSizeInCurrentProcess++;
    }
    size_t sumSizeVectorInPrevProcesses = 0;
    for (size_t i = 0; i < rank; ++i) {
        sumSizeVectorInPrevProcesses += N / sizeProccess;
        if ((N % sizeProccess != 0) && (N % sizeProccess >= i + 1)) {
            sumSizeVectorInPrevProcesses++;
        }
    }
    pieceVectorOfMatrix = calloc(vectorSizeInCurrentProcess * N, sizeof(double));
    for (size_t i = 0; i < vectorSizeInCurrentProcess; ++i) {
        for (size_t j = 0; j < N; ++j) {
            pieceVectorOfMatrix[i * N + j] = 1;
            if (j == i + sumSizeVectorInPrevProcesses) {
                pieceVectorOfMatrix[i * N + j] = 2;
            }
        }
    }

    size_t shiftSize = N / sizeProccess + (N % sizeProccess != 0);

    double *vectorU = calloc(shiftSize, sizeof(double));
    for (size_t i = 0; i < vectorSizeInCurrentProcess; ++i) {
        vectorU[i] = sin(2 * PI * (i + 1 + sumSizeVectorInPrevProcesses) / N);
    }

    double *vectorX = NULL;
    double *vectorB = NULL;
    vectorX = calloc(shiftSize, sizeof(double));
    vectorB = calloc(shiftSize, sizeof(double));

    double *vectorAxn_b = NULL;
    vectorAxn_b = calloc(shiftSize, sizeof(double));

    mulMatrixVector(pieceVectorOfMatrix, vectorU, vectorB,
                    vectorSizeInCurrentProcess, shiftSize, sumSizeVectorInPrevProcesses);

    double startTime = MPI_Wtime();

    double normB = getNorm(vectorB, vectorSizeInCurrentProcess);

    for(size_t k = 0; 1; ++k) {
        setZeroVector(vectorAxn_b, vectorSizeInCurrentProcess);
        mulMatrixVector(pieceVectorOfMatrix, vectorX, vectorAxn_b,
                        vectorSizeInCurrentProcess, shiftSize, sumSizeVectorInPrevProcesses);
        subVector(vectorAxn_b, vectorB, vectorSizeInCurrentProcess);

        double normAx_b = getNorm(vectorAxn_b, vectorSizeInCurrentProcess);
        if (normAx_b / normB < epsilon) {
            break;
        }

        F_ofXn(vectorX, vectorAxn_b, vectorSizeInCurrentProcess);
    }

    double endTime = MPI_Wtime();
    double time = endTime - startTime;

    double finalTime = 0;
    MPI_Reduce(&time, &finalTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("%f\n", finalTime);
    }


    free(pieceVectorOfMatrix);
    free(vectorU);
    free(vectorX);
    free(vectorB);
    free(vectorAxn_b);
    MPI_Finalize();
    return 0;
}