#include "../include/TaskQueue.h"
#include "../include/Worker.h"
#include "../include/Receiver.h"
#include "../include/Sender.h"
#include <mpi.h>
#include <condition_variable>
#include <thread>
#include <iostream>


void printResults(const std::shared_ptr<Worker>& worker, double time, int sumWeight);


int main(int argc, char** argv) {
    bool balancing = true;
    int processID;
    int processCount;
    double beginningTime;
    double endingTime;
    int totalTasksCount = 500;
    int totalSumWeight = 5'000'000;


    int required = MPI_THREAD_MULTIPLE;
    int provided;

    MPI_Init_thread(&argc, &argv, required, &provided);

    if (required != provided) {
        std::cout << "MPI cannot provide required thread support level" << std::endl;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &processID);
    MPI_Comm_size(MPI_COMM_WORLD, &processCount);

    std::shared_ptr<TaskQueue> const taskQueue = std::make_shared<TaskQueue>(totalTasksCount);
    std::shared_ptr<std::mutex> const mutex = std::make_shared<std::mutex>();
    std::shared_ptr<std::condition_variable> const workerCondition = std::make_shared<std::condition_variable>();
    std::shared_ptr<std::condition_variable> const receiverCondition = std::make_shared<std::condition_variable>();


    std::shared_ptr<Worker> worker = std::make_shared<Worker>(processID,
                                                              processCount,
                                                              taskQueue,
                                                              mutex,
                                                              workerCondition,
                                                              receiverCondition,
                                                              totalTasksCount,
                                                              totalSumWeight
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
    if (balancing) {
        std::thread receiverThread (&Receiver::start, receiver);
        std::thread senderThread (&Sender::start, sender);

        worker->start();
        receiverThread.join();
        senderThread.join();
    } else {
        worker->startWithoutBalancing();
    }
    endingTime = MPI_Wtime();
    double tmpTime {endingTime - beginningTime};
    double time {};

    MPI_Reduce(&tmpTime,
               &time,
               1,
               MPI_DOUBLE,
               MPI_MAX,
               0,
               MPI_COMM_WORLD);


    printResults(worker, time, totalSumWeight);

    MPI_Finalize();
    return 0;
}



void printResults(const std::shared_ptr<Worker>& worker, double time, int sumWeight) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (worker->getProcessID() == 0) {
        std::cout << "Total weight of tasks of workers:" << std::endl;
        std::cout << "Worker       | startSumWeight | endSumWeight |" << std::endl;
    }
    for (int i {}; i < worker->getProcessCount(); ++i) {
        if (i == worker->getProcessID()) {
            std::string startSumString = std::to_string(worker->getStartSumWeight());
            std::string endSumString = std::to_string(worker->getEndSumWeight());
            std::string spaces1 = std::string(15 - worker->to_string().size(), ' ');
            std::string spaces2 = std::string(17 - startSumString.size(), ' ');
            std::string spaces3 = std::string(15 - endSumString.size(), ' ');
            std::string line = worker->to_string() + spaces1 + startSumString + spaces2
                    + endSumString;
            std::cout << line << std::endl;
        }
        std::this_thread::sleep_for(std::chrono::nanoseconds(100));

        MPI_Barrier(MPI_COMM_WORLD);
    }

    double tmpStartSumWeight = worker->getStartSumWeight();
    double totalStartSumWeight {};
    double tmpEndSumWeight = worker->getEndSumWeight();
    double totalEndSumWeight {};

    MPI_Allreduce(&tmpStartSumWeight,
               &totalStartSumWeight,
               1,
               MPI_DOUBLE,
               MPI_SUM,
               //0,
               MPI_COMM_WORLD);
    MPI_Allreduce(&tmpEndSumWeight,
               &totalEndSumWeight,
               1,
               MPI_DOUBLE,
               MPI_SUM,
               //0,
               MPI_COMM_WORLD);


    MPI_Barrier(MPI_COMM_WORLD);
    std::this_thread::sleep_for(std::chrono::nanoseconds(1000));
    if (worker->getProcessID() == 0) {
        std::cout << std::endl;
        std::cout << "sumWeight: " << sumWeight << std::endl;
        std::cout << "totalStartSumWeight: " << totalStartSumWeight << std::endl;
        std::cout << "totalEndSumWeight: " << totalEndSumWeight << std::endl;

        if (totalStartSumWeight != sumWeight || totalEndSumWeight != sumWeight) {
            std::cout << "SUMS AREN'T SAME" << std::endl;
            std::cout << "Difference with start: " << sumWeight - totalStartSumWeight << std::endl;
            std::cout << "Difference with end:   " << sumWeight - totalEndSumWeight << std::endl;

        } else {
            std::cout << "sums are same" << std::endl;

        }
        std::cout << "\nTime: " << time << std::endl;
    }
}

