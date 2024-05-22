#include "Sender.h"
#include <mpi/mpi.h>

Sender::Sender(int processID, std::shared_ptr<TaskQueue>& taskQueue, std::shared_ptr<std::mutex>& mutex) {
    this->processID = processID;
    this->taskQueue = taskQueue;
    this->mutex = mutex;
    this->running = true;
}
void Sender::start() {
    while (true) {
        int receiveProcessID;
        Task task = Task();

        printf("Sender %d waiting for request\n", processID);
        MPI_Recv(&receiveProcessID,1,MPI_INT,MPI_ANY_SOURCE,REQUEST_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        if (receiveProcessID == TERMINATION_SIGNAL) {
            printf("Sender %d received termination signal\n", processID);
            break;
        }

        printf("Sender %d received request from process %d\n", processID, receiveProcessID);

        //pthread_mutex_lock(&mutex);
        this->mutex->lock();
        if (!this->taskQueue->isEmpty()) {
            task = this->taskQueue->pop();
            printf("Sender %d sent task %d of process %d to process %d\n", processID, task.getID(), task.getProcessID(), receiveProcessID);
        } else {
            task = Task::createEmptyTask(this->processID);
            printf("Sender %d sent empty queue response to process %d\n", processID, receiveProcessID);
        }
        this->mutex->unlock();
        //pthread_mutex_unlock(&mutex);

        MPI_Send(&task, sizeof(task), MPI_BYTE, receiveProcessID, RESPONSE_TAG, MPI_COMM_WORLD);
    }

    printf("Sender %d finished\n", processID);
}


