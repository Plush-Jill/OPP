#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>


#define N 10
#define EPSILON 0.000001
#define TAU (double)(0.01)
#define MAX_ITERATION_COUNT 1000000000
double vectorBNorm;


void mulMatrixVector(const double *matrix, const double *inputVector, double *outputVector) {
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            outputVector[i] += matrix[i * N + j] * inputVector[j];
        }
    }
}
void subVector(double *vector1, const double *vector2) {
    for (size_t i = 0; i < N; ++i) {
        vector1[i] -= vector2[i];
    }
}
void initMatrixA_v0 (double* matrixA){
    for (int i = 0; i < N; ++i){
        for (int j = 0; j < N; ++j){
            matrixA[i * N + j] = 1;
            if (i == j){
                matrixA[i * N + j] = 2;
            }
        }
    }
}
void initVectorX_v0 (double* vectorX){
    for (int i = 0; i < N; ++i) {
        vectorX[i] = 0;
    }

}
void initVectorB_v0 (double* vectorB){
    for (int i = 0; i < N; ++i) {
        vectorB[i] = N + 1;
    }
}
void calcAxb_v0 (double* vectorAxb, const double* matrixA, const double* vectorX, const double* vectorB){
    /*for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            vectorAxb[i] = matrixA[i * N + j] * vectorX[j];
        }
        vectorAxb[i] -= vectorB[i];
    }*/
    mulMatrixVector(matrixA, vectorX, vectorAxb);
    subVector(vectorAxb, vectorB);
}
void setVectorXtoFofX_v0 (double* vectorX, const double* vectorAxb) {
    for (int i = 0; i < N; ++i){
        vectorX[i] = vectorX[i] - TAU * (vectorAxb[i]);
    }
}
double calcVectorNorm_v0 (const double* vector){
    double res = 0;
    for (size_t i = 0; i < N; ++i) {
        double a = vector[i];
        res += (a * a);
    }
    return res;
}
double calculateSquareRootOfNorm_v0 (const double* vector) {
    return sqrt(calcVectorNorm_v0(vector));
}
void printVector_v0 (double* vector){
    printf("(");
    for (int i = 0; i < N; ++i){
        i == N - 1 ? printf("%lf)\n", vector[i]) : printf("%lf, ", vector[i]);
    }
}
void printMatrix_v0 (double* matrix) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            fprintf(stdout, "%lf ", matrix[i * N + j]);
        }
        fprintf(stdout, "\n");
    }
}
double G_of_Xn (double* matrixA, double* vectorX, double* vectorB, double* vectorAXB) {
    //calcAxb_v0(vectorAXB, matrixA, vectorX, vectorB);

    double upPart = calcVectorNorm_v0(vectorAXB);

    return upPart/vectorBNorm;
}


int main() {
    int iterationsCount;
    double* vectorX = malloc(sizeof(double) * N);
    double* vectorB = malloc(sizeof(double) * N);
    double* vectorAxb = malloc(sizeof(double) * N);
    double* matrixA = malloc(sizeof(double) * N * N);
    double vectorAxbNorm;
    double gOfX;
    double accuracy = N;

    initVectorX_v0(vectorX);
    initVectorB_v0(vectorB);
    initMatrixA_v0(matrixA);
    double beginningTime;
    double endingTime;

    fprintf(stdout, "Matrix:\n");
    printMatrix_v0(matrixA);
    fprintf(stdout, "vectorX:\n");
    printVector_v0(vectorX);
    fprintf(stdout, "vectorB:\n");
    printVector_v0(vectorB);
    fprintf(stdout, "\n\n\n");

    vectorBNorm = calcVectorNorm_v0(vectorB);

    int trash;
    beginningTime = omp_get_wtime();
    for (iterationsCount = 0; iterationsCount < MAX_ITERATION_COUNT; ++iterationsCount){

        /*calcAxb_v0(vectorAxb, matrixA, vectorX, vectorB);
        setVectorXtoFofX_v0(vectorX, vectorAxb);
        gOfX = G_of_Xn(matrixA, vectorX, vectorB, vectorAxb);
        if (gOfX < EPSILON) {
            break;
        }*/
        mulMatrixVector(matrixA, vectorX, vectorAxb);
        subVector(vectorAxb, vectorB);
        vectorAxbNorm = calcVectorNorm_v0(vectorAxb);
        if (vectorAxbNorm / vectorBNorm < EPSILON) {
            // printf("%ld\n", k);
            break;
        }
        for (size_t i = 0; i < N; ++i) {
            vectorX[i] = vectorX[i] - (TAU * vectorAxb[i]);
        }

        if (iterationsCount % 10000 == 0) {
            fprintf(stdout, "vectorX:\n");
            printVector_v0(vectorX);
            fprintf(stderr, "G of Xn: %lf\n", gOfX);
            fprintf(stderr, "Current iter: %d\n", iterationsCount);
        }

    }
    endingTime = omp_get_wtime();
    fprintf(stdout, "Time: %lf sec.\nCompleted iterations: %d\n", endingTime - beginningTime, iterationsCount);


    fprintf(stdout, "vectorX: ");
    printVector_v0(vectorX);

    free(vectorX);
    free(vectorB);
    free(vectorAxb);
    free(matrixA);
}