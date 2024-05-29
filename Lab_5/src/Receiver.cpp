#include "../include/Receiver.h"
#include <mpi.h>
#include <random>
#include <utility>

void Receiver::start() {

    while (this->isRunning()) {
        int receivedTasksCount = 0;
        std::vector<Task> tasks{};
        tasks.reserve(Receiver::maxReceivingTasksCount);

        {
            std::unique_lock<std::mutex> lock(*(this->mutex));
            while (!this->taskQueue->isEmpty()) {
                this->receiverCondition->wait(lock);
            }
        }
        for (int i {processCount - 1}; i >= 0; --i) {
            if (receivedTasksCount > maxReceivingTasksCount) {
                break;
            }
            bool ableToSendMe = isAbleToSendMeTask(i);

            if (!ableToSendMe) {
                continue;
            }
            int taskCountForReceiving {};

            ///sending request for tasks count which current process can get from process i
            MPI_Send(&this->processID,
                     1,
                     MPI_INT,
                     i,
                     Receiver::taskCountRequestMPITag,
                     MPI_COMM_WORLD);
            MPI_Recv(&taskCountForReceiving,
                     1,
                     MPI_INT,
                     i,
                     Receiver::taskCountReplyMPITag,
                     MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
            if (taskCountForReceiving == 0) {
                this->otherProcessesWithTasks.erase(i);
                continue;
            }

            ///receiving taskCountForReceiving tasks from process i
            MPI_Recv(tasks.data(),
                     static_cast<int>(sizeof(Task)) * taskCountForReceiving,
                     MPI_BYTE,
                     i,
                     Receiver::taskReplyMPITag,
                     MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);

            for (int j {}; j < taskCountForReceiving; ++j) {
                    this->mutex->lock();
                    this->taskQueue->push(tasks[j]);
                    this->mutex->unlock();
                    ++receivedTasksCount;
            }
        }

        if (receivedTasksCount == 0) {
            this->stop();
        }

        this->workerCondition->notify_one();
    }

    MPI_Barrier(MPI_COMM_WORLD);
    int endSignal = Receiver::endingSignal;
    MPI_Send(&endSignal,
             1,
             MPI_INT,
             processID,
             Receiver::taskCountRequestMPITag,
             MPI_COMM_WORLD);

}

bool Receiver::isAbleToSendMeTask(int targetProcessID) {
    for (int x : otherProcessesWithTasks) {
        if (x == targetProcessID) {
            return true;
        }
    }
    return false;
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
    for (int i {}; i < processCount; ++i) {
        this->otherProcessesWithTasks.insert(i);
    }
    this->otherProcessesWithTasks.erase(this->processID);

}

std::string Receiver::to_string() const {
    std::string string {};
    string += "[Receiver " + std::to_string(this->processID) + "]";
    return string;
}
