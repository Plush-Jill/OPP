#include "JacobiSolver.h"
#include <cmath>


double JacobiSolver::getX(int i) const {
    return x0 + i * Hx;
}
double JacobiSolver::getY(int i) const {
    return y0 + i * Hy;
}
double JacobiSolver::getZ(int i) const {
    return z0 + i * Hz;
}
int JacobiSolver::getIndex(int x, int y, int z) const {
    return x * Ny * Nz + y * Nz + z;
}
double JacobiSolver::getAccuracyEstimate() const {
    return maxDiff;
}
int JacobiSolver::getProcessRank() const {
    return processRank;
}


double JacobiSolver::calcPhi(double x, double y, double z) const {
    return x * x + y * y + z * z;
}
double JacobiSolver::calcRho(double x, double y, double z) const {
    return 6 - this->A * calcPhi(x, y, z);

}

void JacobiSolver::divideAreaIntoLayers() {
    int offset = 0;

    for (int i = 0; i < this->processCount; ++i) {
        this->layerHeights[i] = this->Nx / this->processCount;
        if (i < this->Nx % this->processCount) {
            this->layerHeights[i]++;
        }
        this->offsets[i] = offset;
        offset += this->layerHeights[i];
    }
}
void JacobiSolver::initLayers() {
    for (int i = 0; i < this->layerHeights[this->processRank]; ++i) {
        for (int j = 0; j < this->Ny; ++j) {
            for (int k = 0; k < this->Nz; ++k) {
                bool isBorder = (this->offsets[this->processRank] + i == 0) || (j == 0) || (k == 0)
                                || (this->offsets[this->processRank] + i == this->Nx - 1)
                                || (j == this->Ny - 1) || (k == this->Nz - 1);

                if (isBorder) {
                    previousFunctionValue[getIndex(i, j, k)] = calcPhi(getX(this->offsets[this->processRank] + i), getY(j), getZ(k));
                    currentFunctionValue[getIndex(i, j, k)] = calcPhi(getX(this->offsets[this->processRank] + i), getY(j), getZ(k));
                } else {
                    previousFunctionValue[getIndex(i, j, k)] = 0;
                    currentFunctionValue[getIndex(i, j, k)] = 0;
                }
            }
        }
    }
}


void JacobiSolver::solveEquation() {
    double maxDifferenceFromCenter;
    double maxDifferenceFromBorder;

    do {
        MPI_Iallreduce(&this->previousProcessMaxDiff, &this->maxDiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, &this->reduceMaxDiffRequest);

        swapPreviousAndCurrentValues();

        sendReceiveBorders();

        maxDifferenceFromCenter = calcCenter();

        waitBorders();

        maxDifferenceFromBorder = calcBorder();

        MPI_Wait(&reduceMaxDiffRequest, MPI_STATUS_IGNORE);
        previousProcessMaxDiff = std::max(maxDifferenceFromCenter, maxDifferenceFromBorder);
    } while (maxDiff >= this->EPSILON);


    swapPreviousAndCurrentValues();

    maxDiff = calcMaxDiff();

}

void JacobiSolver::swapPreviousAndCurrentValues() {
    std::swap(this->previousFunctionValue, this->currentFunctionValue);
}
double JacobiSolver::calcCenter() {
    double f_i;
    double f_j;
    double f_k;
    double tmpMaxDiff;
    double maxDiff_ = 0.0;
    for (int i = 1; i < this->layerHeights[this->processRank] - 1; ++i) {
        for (int j = 1; j < this->Ny - 1; ++j) {
            for (int k = 1; k < this->Nz - 1; ++k) {
                f_i = (this->previousFunctionValue[getIndex(i + 1, j, k)] + this->previousFunctionValue[getIndex(i - 1, j, k)]) / this->Hx2;
                f_j = (this->previousFunctionValue[getIndex(i, j + 1, k)] + this->previousFunctionValue[getIndex(i, j - 1, k)]) / this->Hy2;
                f_k = (this->previousFunctionValue[getIndex(i, j, k + 1)] + this->previousFunctionValue[getIndex(i, j, k - 1)]) / this->Hz2;

                this->currentFunctionValue[getIndex(i, j, k)] =
                        (f_i + f_j + f_k - calcRho(getX(this->offsets[this->processRank] + i), getY(j), getZ(k))) * this->nextIterationConst;

                tmpMaxDiff = std::fabs(this->currentFunctionValue[getIndex(i, j, k)] - this->previousFunctionValue[getIndex(i, j, k)]);
                maxDiff_ = std::max(maxDiff_, tmpMaxDiff);
            }
        }
    }
    return maxDiff_;
}
double JacobiSolver::calcBorder() {
    double f_i;
    double f_j;
    double f_k;
    double tmpMaxDiff;
    double maxDiff_ = 0.0;

    for (int j = 1; j < this->Ny - 1; ++j)
        for (int k = 1; k < this->Nz - 1; ++k) {

            // up border
            if (this->processRank != 0) {
                f_i = (this->previousFunctionValue[getIndex(1, j, k)] + this->upBorderLayer[getIndex(0, j, k)]) / this->Hx2;
                f_j = (this->previousFunctionValue[getIndex(0, j + 1, k)] + this->previousFunctionValue[getIndex(0, j - 1, k)]) / this->Hy2;
                f_k = (this->previousFunctionValue[getIndex(0, j, k + 1)] + this->previousFunctionValue[getIndex(0, j, k - 1)]) / this->Hz2;

                this->currentFunctionValue[getIndex(0, j, k)] = (f_i + f_j + f_k - calcRho(getX(this->offsets[this->processRank]), getY(j), getZ(k))) * this->nextIterationConst;


                tmpMaxDiff = std::abs(this->currentFunctionValue[getIndex(0, j, k)] - previousFunctionValue[getIndex(0, j, k)]);
                maxDiff_ = std::max(maxDiff_, tmpMaxDiff);
            }

            // down border
            if (this->processRank != this->processCount - 1) {
                f_i = (this->previousFunctionValue[getIndex(this->layerHeights[this->processRank] - 2, j, k)] + downBorderLayer[getIndex(0, j, k)]) / this->Hx2;
                f_j = (this->previousFunctionValue[getIndex(this->layerHeights[this->processRank] - 1, j + 1, k)] + previousFunctionValue[getIndex(this->layerHeights[this->processRank] - 1, j - 1, k)]) / this->Hy2;
                f_k = (this->previousFunctionValue[getIndex(this->layerHeights[this->processRank] - 1, j, k + 1)] + previousFunctionValue[getIndex(this->layerHeights[this->processRank] - 1, j, k - 1)]) / this->Hz2;

                this->currentFunctionValue[getIndex(this->layerHeights[this->processRank] - 1, j, k)] = (f_i + f_j + f_k -
                                                                                                         calcRho(getX(this->offsets[this->processRank] + this->layerHeights[this->processRank] - 1), getY(j), getZ(k))) / (2 / this->Hx2 + 2 / this->Hy2 + 2 / this->Hz2 + this->A);

                tmpMaxDiff = std::abs(currentFunctionValue[getIndex(this->layerHeights[this->processRank] - 1, j, k)] - previousFunctionValue[getIndex(this->layerHeights[this->processRank] - 1, j, k)]);
                maxDiff_ = std::max(maxDiff_, tmpMaxDiff);
            }
        }

    return maxDiff_;
}
double JacobiSolver::calcMaxDiff() {
    double module;
    double maxProcessDelta = 0.0;
    double maxDelta = 0.0;

    for (int i = 0; i < this->layerHeights[this->processRank]; ++i) {
        for (int j = 0; j < this->Ny; ++j) {
            for (int k = 0; k < this->Nz; ++k) {
                module = fabs(this->currentFunctionValue[getIndex(i, j, k)] - calcPhi(getX(this->offsets[this->processRank] + i), getY(j), getZ(k)));
                maxProcessDelta = std::max(maxProcessDelta, module);
            }
        }
    }

    MPI_Allreduce(&maxProcessDelta, &maxDelta, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return maxDelta;
}
JacobiSolver::JacobiSolver() : x0(0), y0(0), z0(0),
                               Dx(2.0), Dy(2.0), Dz(2.0),
                               Nx(200), Ny(200), Nz(200),
                               Hx(Dx / (Nx - 1)), Hy(Dy / (Ny - 1)), Hz(Dz / (Nz - 1)),
                               Hx2(Hx * Hx), Hy2(Hy * Hy), Hz2(Hz * Hz),
                               A(1.0E5), EPSILON(1.0E-4), nextIterationConst(1 / ((2 / this->Hx2) + (2 / this->Hy2) + (2 / this->Hz2) + A))
                               {
    MPI_Comm_size(MPI_COMM_WORLD, &this->processCount);
    MPI_Comm_rank(MPI_COMM_WORLD, &this->processRank);

    this->maxDiff = 0;
    this->previousProcessMaxDiff = this->EPSILON;

    this->layerHeights = std::vector<int>(this->processCount);
    this->offsets = std::vector<int>(this->processCount);

    divideAreaIntoLayers();

    this->previousFunctionValue = std::vector<double>(this->layerHeights[this->processRank] * Ny * Nz);
    this->currentFunctionValue = std::vector<double>(this->layerHeights[this->processRank] * Ny * Nz);

    initLayers();

    this->upBorderLayer = std::vector<double>(Ny * Nz);
    this->downBorderLayer = std::vector<double>(Ny * Nz);

}

void JacobiSolver::sendReceiveBorders() {
    if (this->processRank != 0) {
        MPI_Isend(this->previousFunctionValue.data(), this->Ny * this->Nz, MPI_DOUBLE, this->processRank - 1, this->processRank, MPI_COMM_WORLD, &this->sendUpRequest);
        MPI_Irecv(this->upBorderLayer.data(), this->Ny * this->Nz, MPI_DOUBLE, this->processRank - 1, this->processRank - 1, MPI_COMM_WORLD, &this->receiveUpRequest);
    }

    if (this->processRank != this->processCount - 1) {
        double* previousDownBorder = this->previousFunctionValue.data() + (this->layerHeights[processRank] - 1) * this->Ny * this->Nz;
        MPI_Isend(previousDownBorder, this->Ny * this->Nz, MPI_DOUBLE, this->processRank + 1, this->processRank, MPI_COMM_WORLD, &this->sendDownRequest);
        MPI_Irecv(this->downBorderLayer.data(), this->Ny * this->Nz, MPI_DOUBLE, this->processRank + 1, this->processRank + 1, MPI_COMM_WORLD, &this->receiveDownRequest);
    }
}

void JacobiSolver::waitBorders() {
    if (processRank != 0) {
        MPI_Wait(&sendUpRequest, MPI_STATUS_IGNORE);
        MPI_Wait(&receiveUpRequest, MPI_STATUS_IGNORE);
    }

    if (processRank != processCount - 1) {
        MPI_Wait(&sendDownRequest, MPI_STATUS_IGNORE);
        MPI_Wait(&receiveDownRequest, MPI_STATUS_IGNORE);
    }
}
