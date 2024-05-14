#include <vector>
#include <mpi.h>

class JacobiSolver {
private:
    const double x0;
    const double y0;
    const double z0;

    const double Dx;
    const double Dy;
    const double Dz;

    const int Nx;
    const int Ny;
    const int Nz;

    const double Hx;
    const double Hy;
    const double Hz;

    const double Hx2;
    const double Hy2;
    const double Hz2;

    const double A;
    const double EPSILON;
    const double nextIterationConst;

    int processRank;
    int processCount;

    std::vector<int> layerHeights;
    std::vector<int> offsets;

    std::vector<double> upBorderLayer;
    std::vector<double> downBorderLayer;

    std::vector<double> previousFunctionValue;
    std::vector<double> currentFunctionValue;

    double previousProcessMaxDiff;
    double maxDiff;

    MPI_Request sendUpRequest;
    MPI_Request sendDownRequest;
    MPI_Request receiveUpRequest;
    MPI_Request receiveDownRequest;
    MPI_Request reduceMaxDiffRequest;

    [[nodiscard]] double getX(int i) const;
    [[nodiscard]] double getY(int j) const;
    [[nodiscard]] double getZ(int k) const;
    [[nodiscard]] int getIndex(int x, int y, int z) const;

    void divideAreaIntoLayers();
    void initLayers();

    [[nodiscard]] double calcPhi(double x, double y, double z) const;
    [[nodiscard]] double calcRho(double x, double y, double z) const;

    void swapPreviousAndCurrentValues();
    double calcCenter();
    double calcBorder();
    double calcMaxDiff();

public:
    JacobiSolver();
    [[nodiscard]] double getAccuracyEstimate() const;
    [[nodiscard]] int getProcessRank() const;
    void solveEquation();
};


