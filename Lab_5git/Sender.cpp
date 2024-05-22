#include "Sender.h"
#include <mpi.h>
#include <thread>
#include <utility>

void Sender::start() {
    while (true) {
        int receiverProcessID;
        Task task {};

        MPI_Recv(&receiverProcessID,
                 1,
                 MPI_INT,
                 MPI_ANY_SOURCE,
                 Sender::taskRequestMPITag,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

        if (receiverProcessID == Sender::endingSignal) {
            break;
        }
        this->mutex->lock();
        if (!this->taskQueue->isEmpty()){
            task = this->taskQueue->pop();
        } else {
            task = Task::createEmptyTask(this->processID);
        }
        this->mutex->unlock();

        MPI_Send(&task,
                 sizeof(task),
                 MPI_BYTE,
                 receiverProcessID,
                 Sender::taskReplyMPITag,
                 MPI_COMM_WORLD);
    }

    std::this_thread::yield();
}

Sender::Sender(int processID,
               int processCount,
               std::shared_ptr<TaskQueue> taskQueue,
               std::shared_ptr<std::mutex> mutex,
               std::shared_ptr<std::condition_variable> workerCondition,
               std::shared_ptr<std::condition_variable> receiverCondition
               ) :
               processID(processID), processCount(processCount),
               taskQueue(std::move(taskQueue)), mutex(std::move(mutex)),
               workerCondition(std::move(workerCondition)), receiverCondition(std::move(receiverCondition)), running(true)
               {

}

void Sender::stop() {
    this->running = false;
}

