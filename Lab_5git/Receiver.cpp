#include "Receiver.h"
#include "Worker.h"
#include "Sender.h"
#include <mpi/mpi.h>

void Receiver::start() {
    while (this->isRunning()) {
        int receivedTasksCount = 0;
        Task task{};
        std::cout << "Receiver " << this->processID << " started waiting" << ", " << "His queue size = " << this->taskQueue->getSize() << std::endl;

        while (!this->taskQueue->isEmpty()) {
            std::unique_lock<std::mutex> lock (*(this->mutex));
            this->receiverCondition->wait(lock);
        }
        std::cout << "Receiver " << this->processID << " ended waiting" << ", " << "His queue size = " << this->taskQueue->getSize() << std::endl;
        this->mutex->unlock();
        this->mutex->lock();
        std::cout << "Receiver " << this->processID << " sends requests" << std::endl;
        this->mutex->unlock();
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
                this->mutex->lock();
                std::cout << "Receiver " << this->processID
                          << " added task " + task.to_string() << std::endl;
                this->mutex->unlock();
                this->mutex->lock();
                this->taskQueue->push(task);
                this->mutex->unlock();

                ++receivedTasksCount;
            }
        }

        if (receivedTasksCount == 0) {

            this->mutex->lock();
            std::cout << "Receiver " << this->processID << " hasn't get valid tasks" << std::endl;
            this->stop();
            this->mutex->unlock();
        }

        this->mutex->lock();
        this->workerCondition->notify_one();
        this->mutex->unlock();
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

Receiver::Receiver(int processID,
                   int processCount,
                   TaskQueue* taskQueue,
                   std::mutex* mutex,
                   std::condition_variable_any* workerCondition,
                   std::condition_variable_any* receiverCondition,
                   Worker* worker,
                   Sender* sender
                   ) :
                   processID(processID), processCount(processCount),
                   taskQueue(taskQueue), mutex(mutex),
                   workerCondition(workerCondition), receiverCondition(receiverCondition),
                   worker(worker), sender(sender),
                   running(true)
                   {

}
