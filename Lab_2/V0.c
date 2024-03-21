#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>


#define N 32
#define EPSILON 0.00001
#define TAU 1E-5
#define MAX_ITERATION_COUNT 100000000


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
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            vectorAxb[i] = matrixA[i * N + j] * vectorX[j];
        }
        vectorAxb[i] -= vectorB[i];
    }
}
void setVectorXtoFofX_v0 (double* vectorX, const double* vectorAxb){
    for (int i = 0; i < N; ++i){
        vectorX[i] -= TAU * (vectorAxb[i]);
    }
}
double calculateEuclideanVectorNorm_v0 (const double* vector){
    double norm = 0;
    for (int i = 0; i < N; ++i){
        norm += vector[i] * vector[i];
    }
    return norm;
}
double calculateSquareRootOfNorm_v0 (const double* vector) {
    return sqrt(calculateEuclideanVectorNorm_v0(vector));
}
void printVector_v0 (double* vector){
    printf("(");
    for (int i = 0; i < N; ++i){
        i == N - 1 ? printf("%lf)\n", vector[i]) : printf("%lf, ", vector[i]);
    }
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

    initVectorX_v0(vectorX);
    initVectorB_v0(vectorB);
    initMatrixA_v0(matrixA);
    double beginningTime;
    double endingTime;

    /*
    fprintf(stdout, "vectorB: ");
    printVector_v0(vectorB);
    fprintf(stdout, "vectorX: ");
    printVector_v0(vectorX);
    fprintf(stdout, "matrixA: ");
    printVector_v0(matrixA);
    */
    //return 0;

    beginningTime = omp_get_wtime();
    for (iterationsCount = 0; accuracy >= EPSILON && iterationsCount < MAX_ITERATION_COUNT; ++iterationsCount){

        calcAxb_v0(vectorAxb, matrixA, vectorX, vectorB);
        //printVector_v0(vectorAxb);
        //return 0;

        setVectorXtoFofX_v0(vectorX, vectorAxb);
        vectorAxbEuclideanNorm = calculateEuclideanVectorNorm_v0(vectorAxb);
        vectorBEuclideanNorm = calculateEuclideanVectorNorm_v0(vectorB);
        accuracy = sqrt(vectorAxbEuclideanNorm / vectorBEuclideanNorm);
        //fprintf(stdout, "Axb norm = %lf, B norm = %lf, accuracy: %lf\n", vectorAxbEuclideanNorm, vectorBEuclideanNorm, accuracy);

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