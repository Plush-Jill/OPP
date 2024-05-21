#include "TaskQueue.h"
#include "Exceptions.h"
#include "Worker.h"
#include "Receiver.h"
#include "Sender.h"
#include <mpi/mpi.h>
#include <condition_variable>
#include <thread>


void initMPI(int &argc, char **&argv);
void printResults(Worker* worker, double beginningTime, double endingTime);


int main(int argc, char** argv) {
    int processID;
    int processCount;
    double beginningTime;
    double endingTime;
    int totalTaskCount = 120;
    int totalSumWeight = 100000;

    try {
        initMPI(argc, argv);
    }catch (InitMPIException& exception) {
        std::cerr << exception.what() << std::endl;
        return 1;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &processID);
    MPI_Comm_size(MPI_COMM_WORLD, &processCount);

    TaskQueue* taskQueue = new TaskQueue(totalTaskCount);
    std::mutex* mutex = new std::mutex();
    std::condition_variable_any* workerCondition = new std::condition_variable_any();
    std::condition_variable_any* receiverCondition = new std::condition_variable_any();


    Worker worker = Worker(processID,
                           processCount,
                           taskQueue,
                           mutex,
                           workerCondition,
                           receiverCondition,
                           totalTaskCount,
                           totalSumWeight
                           );
    Sender sender = Sender(processID,
                           processCount,
                           taskQueue,
                           mutex,
                           workerCondition,
                           receiverCondition
                           );
    Receiver receiver = Receiver(processID,
                                 processCount,
                                 taskQueue,
                                 mutex,
                                 workerCondition,
                                 receiverCondition,
                                 &worker,
                                 &sender
    );
    MPI_Barrier(MPI_COMM_WORLD);
    beginningTime = MPI_Wtime();
    std::thread workerThread (&Worker::start, &worker);
    std::thread receiverThread (&Receiver::start, &receiver);
    std::thread senderThread (&Sender::start, &sender);


    workerThread.join();
    receiverThread.join();
    senderThread.join();

    endingTime = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    if (processID == 0) {
        std::cout << "END" << std::endl;
    }
    printResults(&worker, beginningTime, endingTime);

    MPI_Finalize();
    return 0;
}



void printResults(Worker* worker, double beginningTime, double endingTime) {
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i {}; i < worker->getProcessCount(); ++i) {
        if (i == worker->getProcessID()) {
            std::cout << "Summary weight of worker " << worker->getProcessID() << ":" << std::endl; ///!!!!!!!!!!!!!
            std::cout << "At the beginning: " << worker->getStartSumWeight() << std::endl;
            std::cout << "At the ending: " << worker->getEndSumWeight() << std::endl;
        }
        std::this_thread::sleep_for(std::chrono::nanoseconds(100));
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    std::this_thread::sleep_for(std::chrono::nanoseconds(1000));
    if (worker->getProcessID() == 0) {
        std::cout << "\nTime: " << endingTime - beginningTime << std::endl;
    }
}

void initMPI(int &argc, char **&argv) noexcept(false){
    int required = MPI_THREAD_MULTIPLE;
    int provided;

    MPI_Init_thread(&argc, &argv, required, &provided);

    if (required != provided) {
        throw InitMPIException("MPI cannot provide required thread support level");
    }
}
