#include "../include/Sender.h"
#include <mpi.h>
#include <utility>

void Sender::start() {

    while (true) {
        int receiverProcessID;
        std::vector<Task> tasksForSending{};
        tasksForSending.reserve(Sender::maxSendingTasksCount);
        int taskCountForSending {};

        MPI_Recv(&receiverProcessID,
                 1,
                 MPI_INT,
                 MPI_ANY_SOURCE,
                 Sender::taskCountRequestMPITag,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);



        if (receiverProcessID == Sender::endingSignal) {
            break;
        }
        if (isAbleToSendTasks()) {

            this->mutex->lock();
            if (isEnoughTasksToSending()) {
                for (int i{}; i < Sender::maxSendingTasksCount; ++i) {
                    if (!this->taskQueue->isEmpty()) {
                        tasksForSending.push_back(this->taskQueue->pop());
                        ++taskCountForSending;
                    } else {
                        break;
                    }
                }
            } else {
                stopSending();
            }
            this->mutex->unlock();

        }
        MPI_Send(&taskCountForSending,
                 1,
                 MPI_INT,
                 receiverProcessID,
                 Sender::taskCountReplyMPITag,
                 MPI_COMM_WORLD);

        if (taskCountForSending == 0) {
            continue;
        }
        MPI_Send(tasksForSending.data(),
                 static_cast<int>(sizeof(Task)) * taskCountForSending,
                 MPI_BYTE,
                 receiverProcessID,
                 Sender::taskReplyMPITag,
                 MPI_COMM_WORLD);
    }

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
        workerCondition(std::move(workerCondition)), receiverCondition(std::move(receiverCondition)),
        running(true), ableToSendTasks(true)
               {

}

void Sender::stop() {
    this->running = false;
}

std::string Sender::to_string() const {
    std::string string {};
    string += "[Sender " + std::to_string(this->processID) + "]";
    return string;
}

bool Sender::isEnoughTasksToSending() const {
    return this->taskQueue->getRemainsTasksCount() > Sender::limitForPossibilityOfSending;
}

void Sender::stopSending() {
    this->ableToSendTasks = false;
}

bool Sender::isAbleToSendTasks() const {
    return this->ableToSendTasks;
}

