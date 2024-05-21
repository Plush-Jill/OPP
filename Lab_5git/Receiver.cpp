#include "Receiver.h"
#include "Worker.h"
#include "Sender.h"
#include <mpi/mpi.h>

Receiver::Receiver(int processID,
                   int processCount,
                   TaskQueue* taskQueue,
                   std::mutex* mutex,
                   std::condition_variable* workerCondition,
                   std::condition_variable* receiverCondition,
                   Worker* worker,
                   Sender* sender,
                   pthread_mutex_t* mutexC,
                   pthread_cond_t* workerConditionC,
                   pthread_cond_t* receiverConditionC
) : processID(processID), processCount(processCount),
    taskQueue(taskQueue), mutex(mutex),
    workerCondition(workerCondition), receiverCondition(receiverCondition),
    mutexC(mutexC),
    workerConditionC(workerConditionC), receiverConditionC(receiverConditionC)
    {
//    this->processID = processID;
//    this->processCount = processCount;
//    this->taskQueue = taskQueue;
//    this->mutex = mutex;
//    this->workerCondition = workerCondition;
//    this->receiverCondition = receiverCondition;
    this->worker = worker;
    this->sender = sender;
    this->running = true;
}
void Receiver::start() {
    while (this->isRunning()) {
        int receivedTasksCount = 0;
        Task task{};

        pthread_mutex_lock(this->mutexC);
        std::cout << "Receiver " << this->processID << " started waiting" << ", " << "His queue size = " << this->taskQueue->getSize() << std::endl;

        while (!this->taskQueue->isEmpty()) {
            pthread_cond_wait(this->receiverConditionC, this->mutexC);
        }
        std::cout << "Receiver " << this->processID << " ended waiting" << ", " << "His queue size = " << this->taskQueue->getSize() << std::endl;
        pthread_mutex_unlock(this->mutexC);
        pthread_mutex_lock(this->mutexC);
        std::cout << "Receiver " << this->processID << " sends requests" << std::endl;
        pthread_mutex_unlock(this->mutexC);
        for (int i {}; i < processCount; ++i) {
            if (i == this->processID) {
                continue;
            }

            MPI_Send(&processID,
                     1,
                     MPI_INT,
                     i,
                     REQUEST_TAG,
                     MPI_COMM_WORLD);
            MPI_Recv(&task,
                     sizeof(task),
                     MPI_BYTE,
                     i,
                     RESPONSE_TAG,
                     MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);

            if (!task.isEmpty()) {
                pthread_mutex_lock(this->mutexC);
                std::cout << "Receiver " << this->processID
                          << " added task " + task.to_string() << std::endl;
                pthread_mutex_unlock(this->mutexC);
                pthread_mutex_lock(this->mutexC);
                this->taskQueue->push(task);
                pthread_mutex_unlock(this->mutexC);

                ++receivedTasksCount;
            }
        }

        if (receivedTasksCount == 0) {

            pthread_mutex_lock(this->mutexC);
            std::cout << "Receiver " << this->processID << " hasn't get valid tasks" << std::endl;
            this->stop();
            pthread_mutex_unlock(this->mutexC);
        }

        pthread_mutex_lock(this->mutexC);
        pthread_cond_signal(this->workerConditionC);
        pthread_mutex_unlock(this->mutexC);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    int terminationSignal = TERMINATION_SIGNAL;
    MPI_Send(&terminationSignal,
             1,
             MPI_INT,
             processID,
             REQUEST_TAG,
             MPI_COMM_WORLD);

    std::this_thread::yield();
}

bool Receiver::isRunning() const {
    return this->running;
}

void Receiver::stop() {
    this->running = false;
    this->worker->stop();
    this->sender->stop();
}
