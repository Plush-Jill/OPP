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
    int totalTaskCount = 4000;

    try {
        initMPI(argc, argv);
    }catch (InitMPIException& exception) {
        std::cerr << exception.what() << std::endl;
        return 1;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &processID);
    MPI_Comm_size(MPI_COMM_WORLD, &processCount);

    std::shared_ptr<TaskQueue> taskQueue = std::make_shared<TaskQueue>(totalTaskCount);
    std::cout << std::to_string(taskQueue->getCapacity()) + " " + std::to_string(taskQueue->getCount())
                + " " + std::to_string(taskQueue->getPopIndex()) << std::endl;

    std::shared_ptr<std::mutex> mutex =
            std::make_shared<std::mutex>();
    std::shared_ptr<std::condition_variable> workerCondition =
            std::shared_ptr<std::condition_variable>();
    std::shared_ptr<std::condition_variable> receiverCondition =
            std::shared_ptr<std::condition_variable>();


    Worker worker = Worker(processID, processCount, totalTaskCount,
                           mutex, workerCondition, receiverCondition, taskQueue);
    Receiver receiver = Receiver(processID, processCount, taskQueue,
                                 mutex, receiverCondition,workerCondition);
    Sender sender = Sender(processID, taskQueue, mutex);


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
