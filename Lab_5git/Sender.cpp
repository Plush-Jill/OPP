#include "Sender.h"
#include <mpi/mpi.h>
#include <thread>

void Sender::start() {
    while (true) {
        int receiveProcessID;
        Task task{};

        MPI_Recv(&receiveProcessID,
                 1,
                 MPI_INT,
                 MPI_ANY_SOURCE,
                 REQUEST_TAG,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

        if (receiveProcessID == TERMINATION_SIGNAL) {
            break;
        }
        pthread_mutex_lock(this->mutexC);
        std::cout << "Sender " << this->processID << " received request for task from process "
                  << receiveProcessID << std::endl;
        pthread_mutex_unlock(this->mutexC);
        pthread_mutex_lock(this->mutexC);
        if (!this->taskQueue->isEmpty()){
            task = this->taskQueue->pop();
        } else {
            task = Task::createEmptyTask(processID);
        }
        pthread_mutex_unlock(this->mutexC);

        MPI_Send(&task,
                 sizeof(task),
                 MPI_BYTE,
                 receiveProcessID,
                 RESPONSE_TAG,
                 MPI_COMM_WORLD);
        pthread_mutex_lock(this->mutexC);
        std::cout << "Sender " << this->processID <<
                     " sent task " + task.to_string() + " to process " << receiveProcessID << std::endl;
        pthread_mutex_unlock(this->mutexC);

    }

    std::this_thread::yield();
}

Sender::Sender(int processID,
               int processCount,
               TaskQueue* taskQueue,
               std::mutex* mutex,
               std::condition_variable* workerCondition,
               std::condition_variable* receiverCondition,
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

    this->running = true;
}

void Sender::stop() {
    this->running = false;
}

