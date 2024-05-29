#include "../include/Worker.h"
#include <utility>
#include <complex>
#include <iostream>

void Worker::start() {

    this->initTasks();

    while (true) {

        executeCurrentTask();
        {

            std::unique_lock<std::mutex> lock (*(this->mutex));
            while (this->taskQueue->isEmpty() && this->isRunning()) {
                this->receiverCondition->notify_one();
                this->workerCondition->wait(lock);
            }
        }

        if (!isRunning()) {
            break;
        }
    }
    executeCurrentTask();
}
void Worker::startWithoutBalancing() {
    this->initTasks();
    executeCurrentTask();
}

void Worker::executeCurrentTask() {
    while (true) {
        Task task{};

        this->mutex->lock();
        if (this->taskQueue->isEmpty()) {
            this->mutex->unlock();
            break;
        }
        if (!isEnoughTasksRemainsForDoOnlyOwnTasks()) {
            this->receiverCondition->notify_one();
        }
        task = this->taskQueue->pop();

        this->mutex->unlock();

        double tmp {};
        for (int i {}; i < task.getWeight(); ++i) {
            for (int j {}; j < Worker::timeScaleToExecuteEachTask; ++j) {
                tmp += sqrt(sqrt(sqrt(sqrt(sqrt(i)))));
            }
        }
        tmp = 1 + 1/tmp;
        this->sumForAvoidingCompilerOptimization += tmp;

        this->endSumWeight += task.getWeight();
    }
}
bool Worker::isRunning() const {
    return this->running;
}
void Worker::initTasks() {
    int minWeight = 2 * this->totalSumWeight / (this->totalTasksCount * (this->processCount + 1));
    int nextTaskID = 1;
    int res = totalSumWeight;

    for (int i {}; i < this->totalTasksCount; ++i){
        int weight = minWeight * (i % this->processCount + 1);
        Task task = Task(nextTaskID, this->processID, weight);
//        if (i % this->processCount == this->processID/*this->processID == this->processCount - 1*/) {
        if (i % (this->processCount / 5 + 1) == this->processID) {
            this->mutex->lock();
            this->taskQueue->push(task);
            this->mutex->unlock();


            ++nextTaskID;
            this->startSumWeight += weight;
        }
        res -= weight;
    }
    if (processID == processCount - 1) {
        Task task = Task(nextTaskID, this->processID, res);
        this->mutex->lock();
        this->taskQueue->push(task);
        this->mutex->unlock();

        ++nextTaskID;
        this->startSumWeight += res;
    }
}

Worker::Worker(int processID,
               int processCount,
               std::shared_ptr<TaskQueue> taskQueue,
               std::shared_ptr<std::mutex> mutex,
               std::shared_ptr<std::condition_variable> workerCondition,
               std::shared_ptr<std::condition_variable> receiverCondition,
               int totalTasksCount,
               int totalSumWeight
               ) :
        processID(processID), processCount(processCount),
        taskQueue(std::move(taskQueue)), mutex(std::move(mutex)),
        workerCondition(std::move(workerCondition)), receiverCondition(std::move(receiverCondition)),
        totalSumWeight(totalSumWeight),
        totalTasksCount(totalTasksCount), running(true),
        startSumWeight(0), endSumWeight(0), sumForAvoidingCompilerOptimization(0),
        thisWorkerStartCount(0)
        {

}
void Worker::stop() {
    this->running = false;
}
int Worker::getStartSumWeight() const {
    return this->startSumWeight;
}
int Worker::getEndSumWeight() const {
    return this->endSumWeight;
}
int Worker::getProcessID() const {
    return this->processID;
}
int Worker::getProcessCount() const {
    return this->processCount;
}

std::string Worker::to_string() const {
    std::string string {};
    string += "[Worker " + std::to_string(this->getProcessID()) + "]";
    return string;
}

bool Worker::isEnoughTasksRemainsForDoOnlyOwnTasks() const {
    return this->taskQueue->getRemainsTasksCount() > Worker::taskCountLimitBeforeNotifyingReceiver;
}

