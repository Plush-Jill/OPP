#include "Receiver.h"
#include <mpi/mpi.h>

Receiver::Receiver(int processID, int processCount, std::shared_ptr<TaskQueue> &taskQueue,
                   std::shared_ptr<std::mutex> &mutex,
                   std::shared_ptr<std::condition_variable> &receiverCondition,
                   std::shared_ptr<std::condition_variable> &workerCondition) {
    this->processID = processID;
    this->taskQueue = taskQueue;
    this->mutex = mutex;
    this->workerCondition = workerCondition;
    this->receiverCondition = receiverCondition;
    this->processCount = processCount;
    this->running = true;
}

void Receiver::start() {
    while (isRunning()) {
        int receivedTasks = 0;
        Task task = Task();

        this->mutex->lock();
        if (!this->taskQueue->isEmpty()) {
            std::unique_lock<std::mutex> lock (*(this->mutex));
            receiverCondition->wait(lock);
//            pthread_cond_wait(&receiverCondition, &mutex);
        }
        this->mutex->unlock();

        for (int i = 0; i < this->processCount; ++i) {
            if (i == this->processID) {
                continue;
            }

            printf("Receiver %d sent request to process %d\n", processID, i);
            MPI_Send(&processID, 1, MPI_INT, i, REQUEST_TAG, MPI_COMM_WORLD);
            MPI_Recv(&task, sizeof(task), MPI_BYTE, i, RESPONSE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            if (!task.isEmpty()) {
                printf("Receiver %d received task %d from process %d\n", processID, task.getID(), i);

                this->mutex->lock();
                this->taskQueue->push(task);
                this->mutex->unlock();

                ++receivedTasks;
            } else {
                printf("Receiver %d received empty queue response from process %d\n", processID, i);
            }
        }

        if (receivedTasks == 0) {
            this->mutex->lock();
            this->running = false;
            this->mutex->unlock();
        }

        this->mutex->lock();
        workerCondition->notify_one();
        //pthread_cond_signal(&workerCondition);
        this->mutex->unlock();
    }
}

bool Receiver::isRunning() const {
    return this->running;
}
