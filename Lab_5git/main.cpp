#include "TaskQueue.h"
#include "Exceptions.h"
#include "Worker.h"
#include "Receiver.h"
#include "Sender.h"
#include <mpi/mpi.h>
#include <condition_variable>
#include <thread>


void initMPI(int &argc, char **&argv);
void printResults(int processID, double beginningTime, double endingTime);


int main(int argc, char** argv) {
    int processID;
    int processCount;
    double beginningTime;
    double endingTime;
    int totalTaskCount = 21;
    int totalSumWeight = 1200;

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
    std::condition_variable* workerCondition = new std::condition_variable();
    std::condition_variable* receiverCondition = new std::condition_variable();

    pthread_mutex_t mutexC;
    pthread_cond_t workerConditionC;
    pthread_cond_t receiverConditionC;

    pthread_mutex_init(&mutexC, nullptr);
    pthread_cond_init(&workerConditionC, nullptr);
    pthread_cond_init(&receiverConditionC, nullptr);


    Worker worker = Worker(processID,
                           processCount,
                           taskQueue,
                           mutex,
                           workerCondition,
                           receiverCondition,
                           totalTaskCount,
                           totalSumWeight,
                           &mutexC,
                           &workerConditionC,
                           &receiverConditionC
                           );
    Sender sender = Sender(processID,
                           processCount,
                           taskQueue,
                           mutex,
                           workerCondition,
                           receiverCondition,
                           &mutexC,
                           &workerConditionC,
                           &receiverConditionC
                           );
    Receiver receiver = Receiver(processID,
                                 processCount,
                                 taskQueue,
                                 mutex,
                                 workerCondition,
                                 receiverCondition,
                                 &worker,
                                 &sender,
                                 &mutexC,
                                 &workerConditionC,
                                 &receiverConditionC
    );
    if (processID == 0) {
        std::cout << "Starting..." << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    beginningTime = MPI_Wtime();
    std::thread workerThread (&Worker::start, &worker);
    std::thread receiverThread (&Receiver::start, &receiver);
    std::thread senderThread (&Sender::start, &sender);


    workerThread.join();
    receiverThread.join();
    senderThread.join();

    endingTime = MPI_Wtime();

    printResults(processID, beginningTime, endingTime);

    MPI_Finalize();
    return 0;
}



void printResults(int processID, double beginningTime, double endingTime) {
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "Summary weight " << processID << ": " << std::endl; ///!!!!!!!!!!!!!
    MPI_Barrier(MPI_COMM_WORLD);
    if (processID == 0) {
        std::cout << "Time: " << endingTime - beginningTime << std::endl;
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
