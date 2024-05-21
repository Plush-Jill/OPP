#include "Worker.h"
#include <mpi/mpi.h>
#include <complex>

void Worker::start() {
    pthread_mutex_lock(this->mutexC);
    std::cout << "Worker " << this->processID << " inits tasks." << std::endl;
    pthread_mutex_unlock(this->mutexC);
    this->initTasks();
    pthread_mutex_lock(this->mutexC);
    std::cout << "Worker " << this->processID << " has "
              << this->taskQueue->getSize() << " tasks." << std::endl;
    pthread_mutex_unlock(this->mutexC);
    //MPI_Barrier(MPI_COMM_WORLD);

    while (true) {

        executeCurrentTask();
        pthread_mutex_lock(this->mutexC);
        std::cout << "Worker " << this->processID << " finished current tasks" << std::endl;
        pthread_mutex_unlock(this->mutexC);
        pthread_mutex_lock(this->mutexC);
        while (this->taskQueue->isEmpty() && this->isRunning()) {
            ///receiverCondition->notify_one();
            std::cout << "Worker " << processID << " waiting for tasks" << std::endl;
            ///std::unique_lock<std::mutex> lock(*this->mutex);

            std::cout << "Worker " << processID << " created his std::unique_lock" << std::endl;
            ///workerCondition->wait(lock);
            pthread_cond_signal(this->receiverConditionC);
            pthread_cond_wait(this->workerConditionC, this->mutexC);
        }

        if (!isRunning()) {
            pthread_mutex_unlock(this->mutexC);
            break;
        }
        pthread_mutex_unlock(this->mutexC);
    }
    pthread_mutex_lock(this->mutexC);
    std::cout << "Worker " << this->processID << "finished" << std::endl;
    pthread_mutex_unlock(this->mutexC);
    std::this_thread::yield();
}

void Worker::executeCurrentTask() {
    while (true) {
        Task task{};

        pthread_mutex_lock(this->mutexC);
        if (this->taskQueue->isEmpty()) {
            pthread_mutex_unlock(this->mutexC);
            break;
        }
        task = this->taskQueue->pop();
        pthread_mutex_unlock(this->mutexC);
        pthread_mutex_lock(this->mutexC);
        std::cout << "Worker " << this->processID << " doing task " + task.to_string()
                << ", His queue size was = " << this->taskQueue->getSize() + 1 << std::endl;
        pthread_mutex_unlock(this->mutexC);
        //std::this_thread::sleep_for(std::chrono::nanoseconds(task.getWeight()));
        double res {};
        for (int i = 0; i < task.getWeight(); ++i) {
            for (int j = 0; j < 2000; ++j) {
                res += sqrt(sqrt(sqrt(std::sqrt(i))));
            }
        }

        pthread_mutex_lock(this->mutexC);
        std::cout << "Worker " << this->processID << " did task " + task.to_string()
        << ", His queue size = " << this->taskQueue->getSize() << std::endl;
        pthread_mutex_unlock(this->mutexC);
        this->currentProcessEndSumWeight += task.getWeight();
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
            this->currentProcessStartSumWeight += weight;
        }
    }
}

Worker::Worker(int processID,
               int processCount,
               TaskQueue* taskQueue,
               std::mutex* mutex,
               std::condition_variable* workerCondition,
               std::condition_variable* receiverCondition,
               int taskCount,
               int totalSumWeight,
               pthread_mutex_t* mutexC,
               pthread_cond_t* workerConditionC,
               pthread_cond_t* receiverConditionC
               ) : processID(processID), processCount(processCount),
                    taskQueue(taskQueue), mutex(mutex),
                    workerCondition(workerCondition), receiverCondition(receiverCondition),
                    totalSumWeight(totalSumWeight),
                    mutexC(mutexC),
                    workerConditionC(workerConditionC), receiverConditionC(receiverConditionC)
                {
    this->taskCount = taskCount;
    this->running = true;
    this->currentProcessStartSumWeight = 0;
    this->currentProcessEndSumWeight = 0;
}

void Worker::stop() {
    this->running = false;
}
