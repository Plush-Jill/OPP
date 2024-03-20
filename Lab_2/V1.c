#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define N 15000
#define EPSILON 1E-7
#define TAU 1E-5
#define MAX_ITERATION_COUNT 100000


void initMatrixA_v1 (double* matrixA){
    for (int i = 0; i < N; ++i){
        for (int j = 0; j < N; ++j){
            matrixA[i * N + j] = 1;
            if (i == j){
                matrixA[i * N + j] = 2;
            }
        }
    }
}
void initVectorX_v1 (double* vectorX){
    memset(vectorX, 0, N * sizeof(double));
}
void initVectorB1_v1 (double* vectorB){
    memset(vectorB, N + 1, N * sizeof(double));
}
void calcAxb_v1 (double* vectorAxb, const double* matrixA, const double* vectorX, const double* vectorB){
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            vectorAxb[i] += matrixA[i * N + j] * vectorX[j];
        }
        vectorAxb[i] -= vectorB[i];
    }
}
void setVectorXtoFofX_v1 (double* vectorX, const double* vectorAxb){
    for (int i = 0; i < N; ++i){
        vectorX[i] -= TAU * (vectorAxb[i]);
    }
}
double calculateEuclideanVectorNorm_v1 (const double* vector){
    double norm = 0;
    for (int i = 0; i < N; ++i){
        norm += vector[i] * vector[i];
    }
    return norm;
}
double calculateSquareRootOfNorm_v1 (const double* vector) {
    return sqrt(calculateEuclideanVectorNorm_v1(vector));
}
void printVector_v1 (double* vector){
    printf("(");
    for (int i = 0; i < N; ++i){
        i == N - 1 ? printf("%lf)", vector[i]) : printf("%lf, ", vector[i]);
    }
    //printf(")\n");
}


int main_v1() {
    int iterationsCount;
    double* vectorX = malloc(sizeof(double) * N);
    double* vectorB = malloc(sizeof(double) * N);
    double* vectorAxb = malloc(sizeof(double) * N);
    double* matrixA = malloc(sizeof(double) * N * N);
    double vectorAxbEuclideanNorm;
    double vectorBEuclideanNorm;
    double accuracy = N;

    initVectorX_v1(vectorX);
    initVectorB1_v1(vectorB);
    initMatrixA_v1(matrixA);
    double beginningTime;
    double endingTime;


    for (iterationsCount = 0; accuracy >= EPSILON && iterationsCount < MAX_ITERATION_COUNT; ++iterationsCount){

        calcAxb_v1(vectorAxb, matrixA, vectorX, vectorB);
        setVectorXtoFofX_v1(vectorX, vectorAxb);
        vectorAxbEuclideanNorm = calculateEuclideanVectorNorm_v1(vectorAxb);
        vectorBEuclideanNorm = calculateEuclideanVectorNorm_v1(vectorB);
        accuracy = sqrt(vectorAxbEuclideanNorm / vectorBEuclideanNorm);

    }
    printVector_v1(vectorX);

    free(vectorX);
    free(vectorB);
    free(vectorAxb);
    free(matrixA);
}