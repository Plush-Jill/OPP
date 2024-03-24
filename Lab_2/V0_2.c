#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <omp.h>

#define N 10
#define EPSILON 0.000001
#define TAU (double)(0.01)
#define MAX_ITERATIONS_COUNT 1000000000

void mulMatrixVector_v02(const double* matrix, const double* vector, double* resultVector) {
    for (int i = 0; i < N; ++i) {
        resultVector[i] = 0;
        for (int j = 0; j < N; ++j) {
            resultVector[i] += matrix[i * N + j] * vector[j];
        }
    }
}
void subVector_v02(double* vector1, double* vector2, double* resultVector) {
    for (int i = 0; i < N; ++i) {
        resultVector[i] = vector1[i] - vector2[i];
    }
}
double getNorm(const double *vector) {
    double res = 0;
    for (size_t i = 0; i < N; ++i) {
        double a = vector[i];
        res += (a * a);
    }
    return sqrt(res);
}
void setXtoFofX (double* vectorX, const double* vectorAXB) {
    for (int i = 0; i < N; ++i) {
        vectorX[i] -= TAU * vectorAXB[i];
    }
}
void initMatrixA_v02 (double* matrixA) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            matrixA[i * N + j] = 1;
            if (i == j) {
                matrixA[i * N + j] = 2;
            }
        }
    }
}
void initVectorX_v02 (double* vectorX) {
    for (int i = 0; i < N; ++i) {
        vectorX[i] = 0;
    }

}
void initVectorB_v02 (double* vectorB) {
    for (int i = 0; i < N; ++i) {
        vectorB[i] = N + 1;
    }
}
void printVector_v02 (double* vector) {
    printf("(");
    for (int i = 0; i < N; ++i) {
        i == N - 1 ? printf("%lf)\n", vector[i]) : printf("%lf, ", vector[i]);
    }
}
bool checkResultVector (double* vectorX) {
    double preEpsilon = 1 - EPSILON;
    double postEpsilon = 1 + EPSILON;

    for (int i = 0; i < N; ++i) {
        if (!(vectorX[i] < postEpsilon && vectorX[i] > preEpsilon)) {
            return false;
        }
    }

    return true;
}

int main() {
    double* matrixA = malloc(sizeof(double) * N * N);
    double* vectorB = malloc(N * sizeof(double));
    double* vectorX = malloc(N * sizeof(double));
    double* vectorAX = malloc(N * sizeof(double));
    double* vectorAXB = malloc(N * sizeof(double));

    initMatrixA_v02(matrixA);
    initVectorX_v02(vectorX);
    initVectorB_v02(vectorB);
    double vectorBNorm = getNorm(vectorB);


    double beginningTime = omp_get_wtime();



    int i = 0;
    double vectorAXBNorm;
    double accuracy;
    for (; i <= MAX_ITERATIONS_COUNT; ++i) {
        mulMatrixVector_v02(matrixA, vectorX, vectorAX);
        subVector_v02(vectorAX, vectorB, vectorAXB);
        setXtoFofX(vectorX, vectorAXB);

        vectorAXBNorm = getNorm(vectorAXB);
        accuracy = vectorAXBNorm / vectorBNorm;
        if (accuracy < EPSILON) {
            break;
        }
    }



    double endingTime = omp_get_wtime();
    fprintf(stdout, "Result accuracy: %lf\n", accuracy);
    fprintf(stdout, "Time = %.10lf sec\n", endingTime - beginningTime);
    fprintf(stdout, "Total iterations count: %d\n", i);


    //printVector_v02(vectorX);
    /*if (checkResultVector(vectorX)) {
        printf("Vector is calculated correctly.\n");
    } else {
        printf("Vector is calculated incorrectly.\n");
    }*/



    free(matrixA);
    free(vectorB);
    free(vectorX);
    free(vectorAX);
    free(vectorAXB);

    return 0;
}
