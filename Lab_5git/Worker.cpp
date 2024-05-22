#include "Worker.h"
#include <utility>
#include <complex>
#include <thread>

void Worker::start() {
    this->initTasks();

    while (true) {

        executeCurrentTask();
        while (this->taskQueue->isEmpty() && this->isRunning()) {
            std::unique_lock<std::mutex> lock (*(this->mutex));
            this->receiverCondition->notify_one();
            this->workerCondition->wait(lock);
        }

        if (!isRunning()) {
            this->mutex->unlock();
            break;
        }
        this->mutex->unlock();
    }

    std::this_thread::yield();
}

void Worker::executeCurrentTask() {
    while (true) {
        Task task{};

        this->mutex->lock();
        if (this->taskQueue->isEmpty()) {
            this->mutex->unlock();
            break;
        }
        task = this->taskQueue->pop();
        this->mutex->unlock();

        std::this_thread::sleep_for(std::chrono::nanoseconds(task.getWeight() * 10000));
        this->endSumWeight += task.getWeight();
    }
}
bool Worker::isRunning() const {
    return this->running;
}
void Worker::initTasks() {
    int minWeight = 2 * this->totalSumWeight / (this->taskCount * (this->processCount + 1));
    int nextTaskID = 1;

    for (int i {}; i < this->taskCount; ++i){
        int weight = minWeight * (i % this->processCount + 1);
        Task task = Task(nextTaskID, this->processID, weight);

        if (i % this->processCount == this->processID) {
            this->taskQueue->push(task);
            ++nextTaskID;
            this->startSumWeight += weight;
        }
    }
}

Worker::Worker(int processID,
               int processCount,
               std::shared_ptr<TaskQueue> taskQueue,
               std::shared_ptr<std::mutex> mutex,
               std::shared_ptr<std::condition_variable> workerCondition,
               std::shared_ptr<std::condition_variable> receiverCondition,
               int taskCount,
               int totalSumWeight
               ) :
               processID(processID), processCount(processCount),
               taskQueue(std::move(taskQueue)), mutex(std::move(mutex)),
               workerCondition(std::move(workerCondition)), receiverCondition(std::move(receiverCondition)),
               totalSumWeight(totalSumWeight),
               taskCount(taskCount), running(true),
               startSumWeight(0), endSumWeight(0)
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
