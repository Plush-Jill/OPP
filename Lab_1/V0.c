#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define N 15000
#define EPSILON 1E-7
#define TAU 1E-5
#define MAX_ITERATION_COUNT 100000


void initMatrixA (double* matrixA){
    for (int i = 0; i < N; ++i){
        for (int j = 0; j < N; ++j){
            matrixA[i * N + j] = 1;
            if (i == j){
                matrixA[i * N + j] = 2;
            }
        }
    }
}
void initVectorX (double* vectorX){
    memset(vectorX, 0, N * sizeof(double));
}
void initVectorB (double* vectorB){
    memset(vectorB, N + 1, N * sizeof(double));
}
void calcAxb (double* vectorAxb, const double* matrixA, const double* vectorX, const double* vectorB){
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            vectorAxb[i] += matrixA[i * N + j] * vectorX[j];
        }
        vectorAxb[i] -= vectorB[i];
    }
}
void setVectorXtoFofX (double* vectorX, const double* vectorAxb){
    for (int i = 0; i < N; ++i){
        vectorX[i] -= TAU * (vectorAxb[i]);
    }
}
double calculateEuclideanVectorNorm (const double* vector){
    double norm = 0;
    for (int i = 0; i < N; ++i){
        norm += vector[i] * vector[i];
    }
    return norm;
}
double calculateSquareRootOfNorm (const double* vector) {
    return sqrt(calculateEuclideanVectorNorm(vector));
}
void printVector (double* vector){
    printf("(");
    for (int i = 0; i < N; ++i){
        i == N - 1 ? printf("%lf)", vector[i]) : printf("%lf, ", vector[i]);
    }
    //printf(")\n");
}


int main() {
    int iterationsCount;
    double* vectorX = malloc(sizeof(double) * N);
    double* vectorB = malloc(sizeof(double) * N);
    double* vectorAxb = malloc(sizeof(double) * N);
    double* matrixA = malloc(sizeof(double) * N * N);
    double vectorAxbEuclideanNorm;
    double vectorBEuclideanNorm;
    double accuracy = N;

    initVectorX(vectorX);
    initVectorB(vectorB);
    initMatrixA(matrixA);
    double beginningTime;
    double endingTime;


    for (iterationsCount = 0; accuracy >= EPSILON && iterationsCount < MAX_ITERATION_COUNT; ++iterationsCount){

        calcAxb(vectorAxb, matrixA, vectorX, vectorB);
        setVectorXtoFofX(vectorX, vectorAxb);
        vectorAxbEuclideanNorm = calculateEuclideanVectorNorm(vectorAxb);
        vectorBEuclideanNorm = calculateEuclideanVectorNorm(vectorB);
        accuracy = sqrt(vectorAxbEuclideanNorm / vectorBEuclideanNorm);

    }
    printVector(vectorX);

    free(vectorX);
    free(vectorB);
    free(vectorAxb);
    free(matrixA);
}