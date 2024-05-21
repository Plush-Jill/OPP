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
        this->mutex->lock();
        std::cout << "Sender " << this->processID << " received request for task from process "
                  << receiveProcessID << std::endl;
        this->mutex->unlock();
        this->mutex->lock();
        if (!this->taskQueue->isEmpty()){
            task = this->taskQueue->pop();
        } else {
            task = Task::createEmptyTask(processID);
        }
        this->mutex->unlock();

        MPI_Send(&task,
                 sizeof(task),
                 MPI_BYTE,
                 receiveProcessID,
                 RESPONSE_TAG,
                 MPI_COMM_WORLD);
        this->mutex->lock();
        std::cout << "Sender " << this->processID <<
                     " sent task " + task.to_string() + " to process " << receiveProcessID << std::endl;
        this->mutex->unlock();

    }

    std::this_thread::yield();
}

Sender::Sender(int processID,
               int processCount,
               TaskQueue* taskQueue,
               std::mutex* mutex,
               std::condition_variable_any* workerCondition,
               std::condition_variable_any* receiverCondition
               ) :
               processID(processID), processCount(processCount),
               taskQueue(taskQueue), mutex(mutex),
               workerCondition(workerCondition), receiverCondition(receiverCondition), running(true)
               {

}

void Sender::stop() {
    this->running = false;
}

