#include <iostream>
#include <mpi.h>
#include "JacobiSolver.h"

int main(int argc, char** argv) {
    double beginTime;
    double endTime;
    MPI_Init(&argc, &argv);


    JacobiSolver* jacobiSolver = new JacobiSolver();



    beginTime = MPI_Wtime();

    jacobiSolver->solveEquation();

    endTime = MPI_Wtime();



    if (jacobiSolver->getProcessRank() == 0) {
        printf("Time: %lf\n", endTime - beginTime);
        printf("Max Difference: %lf\n", jacobiSolver->getAccuracyEstimate());
    }



    delete jacobiSolver;
    MPI_Finalize();
    return 0;
}
