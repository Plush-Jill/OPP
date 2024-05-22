#include "../include/Receiver.h"
#include <mpi.h>
#include <thread>

#include <utility>

void Receiver::start() {
    while (this->isRunning()) {
        int receivedTasksCount = 0;
        Task task{};

        {
            std::unique_lock<std::mutex> lock (*(this->mutex));
            while (!this->taskQueue->isEmpty()) {
                this->receiverCondition->wait(lock);
            }
        }
        //this->mutex->unlock();

        for (int i {this->processCount - 1}; i >= 0; --i) {
            if (i == this->processID) {
                continue;
            }

            MPI_Send(&this->processID,
                     1,
                     MPI_INT,
                     i,
                     Receiver::taskRequestMPITag,
                     MPI_COMM_WORLD);
            MPI_Recv(&task,
                     sizeof(task),
                     MPI_BYTE,
                     i,
                     Receiver::taskReplyMPITag,
                     MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);

            if (!task.isEmpty()) {
                this->mutex->lock();
                this->taskQueue->push(task);
                this->mutex->unlock();

                ++receivedTasksCount;
            }
        }

        if (receivedTasksCount == 0) {

            this->mutex->lock();
            this->stop();
            this->mutex->unlock();
        }

        this->workerCondition->notify_one();
    }

    MPI_Barrier(MPI_COMM_WORLD);
    int buf = Receiver::endingSignal;
    MPI_Send(&buf,
             1,
             MPI_INT,
             processID,
             Receiver::taskRequestMPITag,
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
                   std::shared_ptr<TaskQueue> taskQueue,
                   std::shared_ptr<std::mutex> mutex,
                   std::shared_ptr<std::condition_variable> workerCondition,
                   std::shared_ptr<std::condition_variable> receiverCondition,
                   std::shared_ptr<Worker> worker,
                   std::shared_ptr<Sender> sender
                   ) :
                   processID(processID), processCount(processCount),
                   taskQueue(std::move(taskQueue)), mutex(std::move(mutex)),
                   workerCondition(std::move(workerCondition)), receiverCondition(std::move(receiverCondition)),
                   worker(std::move(worker)), sender(std::move(sender)),
                   running(true)
                   {

}
