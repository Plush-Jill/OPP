#include "TaskQueue.h"
#include "Exceptions.h"
#include "Worker.h"
#include "Receiver.h"
#include "Sender.h"
#include <mpi.h>
#include <condition_variable>
#include <thread>
#include <iostream>


void initMPI(int &argc, char **&argv);
void printResults(const std::shared_ptr<Worker>& worker, double time);


int main(int argc, char** argv) {
    int processID;
    int processCount;
    double beginningTime;
    double endingTime;
    int tasksCount = 2000;
    int sumWeight = 8000000;

    try {
        initMPI(argc, argv);
    }catch (InitMPIException& exception) {
        std::cerr << exception.what() << std::endl;
        return 1;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &processID);
    MPI_Comm_size(MPI_COMM_WORLD, &processCount);

    std::shared_ptr<TaskQueue> const taskQueue = std::make_shared<TaskQueue>(tasksCount);
    std::shared_ptr<std::mutex> const mutex = std::make_shared<std::mutex>();
    std::shared_ptr<std::condition_variable> const workerCondition = std::make_shared<std::condition_variable>();
    std::shared_ptr<std::condition_variable> const receiverCondition = std::make_shared<std::condition_variable>();


    std::shared_ptr<Worker> worker = std::make_shared<Worker>(processID,
                           processCount,
                           taskQueue,
                           mutex,
                           workerCondition,
                           receiverCondition,
                           tasksCount,
                           sumWeight
                           );
    std::shared_ptr<Sender> sender = std::make_shared<Sender>(processID,
                           processCount,
                           taskQueue,
                           mutex,
                           workerCondition,
                           receiverCondition
                           );
    std::shared_ptr<Receiver> receiver = std::make_shared<Receiver>(processID,
                                 processCount,
                                 taskQueue,
                                 mutex,
                                 workerCondition,
                                 receiverCondition,
                                 worker,
                                 sender
    );
    MPI_Barrier(MPI_COMM_WORLD);
    beginningTime = MPI_Wtime();
    std::thread workerThread (&Worker::start, worker);
    std::thread receiverThread (&Receiver::start, receiver);
    std::thread senderThread (&Sender::start, sender);


    workerThread.join();
    receiverThread.join();
    senderThread.join();

    endingTime = MPI_Wtime();
    double tmpTime {endingTime - beginningTime};
    double time {};
    /*MPI_Reduce(&tmpTime,
               &time,
               1,
               MPI_DOUBLE,
               MPI_MAX,
               0,
               MPI_COMM_WORLD);*/


    printResults(worker, time);

    MPI_Finalize();
    return 0;
}



void printResults(const std::shared_ptr<Worker>& worker, double time) {
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i {}; i < worker->getProcessCount(); ++i) {
        if (i == worker->getProcessID()) {
            std::cout << "Summary weight of worker " << worker->getProcessID() << ":\n"
            << "At the beginning: " << worker->getStartSumWeight()
            << "\nAt the ending: " << worker->getEndSumWeight() << std::endl;
        }
        std::this_thread::sleep_for(std::chrono::nanoseconds(100));

        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    std::this_thread::sleep_for(std::chrono::nanoseconds(1000));
    if (worker->getProcessID() == 0) {
        std::cout << "\nTime: " << time << std::endl;
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
