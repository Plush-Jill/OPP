#include "Worker.h"
#include <mpi/mpi.h>
#include <complex>

void Worker::start() {
    this->mutex->lock();
    std::cout << "Worker " << this->processID << " inits tasks." << std::endl;
    this->mutex->unlock();
    this->initTasks();
    this->mutex->lock();
    std::cout << "Worker " << this->processID << " has "
              << this->taskQueue->getSize() << " tasks." << std::endl;
    this->mutex->unlock();

    while (true) {

        executeCurrentTask();
        std::cout << "Worker " << this->processID << " finished current tasks" << std::endl;
        while (this->taskQueue->isEmpty() && this->isRunning()) {
            std::unique_lock<std::mutex> lock (*(this->mutex));
            std::cout << "Worker " << processID << " waiting for tasks" << std::endl;
            this->receiverCondition->notify_one();
            this->workerCondition->wait(lock);

        }

        if (!isRunning()) {
            this->mutex->unlock();
            break;
        }
        this->mutex->unlock();
    }
    this->mutex->lock();
    std::cout << "Worker " << this->processID << "finished" << std::endl;
    this->mutex->unlock();
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
        this->mutex->lock();
        std::cout << "Worker " << this->processID << " doing task " + task.to_string()
                << ", His queue size was = " << this->taskQueue->getSize() + 1 << std::endl;
        this->mutex->unlock();


        std::this_thread::sleep_for(std::chrono::nanoseconds(task.getWeight() * 100000));

        this->mutex->lock();
        std::cout << "Worker " << this->processID << " did task " + task.to_string()
        << ", His queue size = " << this->taskQueue->getSize() << std::endl;
        this->mutex->unlock();
        this->endSumWeight += task.getWeight();
    }
}
bool Worker::isRunning() const {
    return this->running;
}
void Worker::initTasks() {
    int minWeight = 2 * this->totalSumWeight / (this->taskCount * (processCount + 1));
    int nextTaskID = 1;

    for (int i {}; i < this->taskCount; ++i){
        int weight = minWeight * (i % processCount + 1);
        Task task = Task(nextTaskID, processID, weight);

        if (i % this->processCount == this->processID) {
            this->taskQueue->push(task);
            ++nextTaskID;
            this->startSumWeight += weight;
        }
    }
}

Worker::Worker(int processID,
               int processCount,
               TaskQueue* taskQueue,
               std::mutex* mutex,
               std::condition_variable_any* workerCondition,
               std::condition_variable_any* receiverCondition,
               int taskCount,
               int totalSumWeight
               ) :
               processID(processID), processCount(processCount),
               taskQueue(taskQueue), mutex(mutex),
               workerCondition(workerCondition), receiverCondition(receiverCondition),
               totalSumWeight(totalSumWeight),
               taskCount(taskCount), running(true),
               startSumWeight(0), endSumWeight(0)
               {

}
void Worker::stop() {
    this->running = false;
}
int Worker::getStartSumWeight() const {
    return startSumWeight;
}
int Worker::getEndSumWeight() const {
    return endSumWeight;
}
int Worker::getProcessID() const {
    return this->processID;
}
int Worker::getProcessCount() const {
    return this->processCount;
}
