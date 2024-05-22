#include "Worker.h"
#include <mpi/mpi.h>

Worker::Worker(int processID, int processCount, int taskCount, std::shared_ptr<std::mutex> &mutex,
               std::shared_ptr<std::condition_variable> &workerCondition,
               std::shared_ptr<std::condition_variable> &receiverCondition, std::shared_ptr<TaskQueue> &taskQueue) {
    std::cout << "In Worker constructor " << 1 << std::endl;
    this->processID = processID;
    this->taskQueue = taskQueue;
    this->mutex = mutex;
    this->taskCount = taskCount;
    std::cout << "In Worker constructor " << 2 << std::endl;
    this->workerCondition = workerCondition;
    this->receiverCondition = receiverCondition;
    this->running = true;
    std::cout << "In Worker constructor " << 3 << std::endl;
    this->currentProcessSumWeight = 0;
    this->processCount = processCount;
}

void Worker::start() {
    std::cout << "Init tasks..." << std::endl;

    initTasks();
    std::cout << "Tasks inited" << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);

    while (true) {
        std::cout << "Executing tasks..." << std::endl;
        executeCurrentTask();

        //this->mutex->lock();
        std::lock_guard<std::mutex> lockGuard(*(this->mutex));
        if (this->taskQueue->isEmpty() && isRunning()){
//            std::lock_guard<std::mutex> lock1(*(this->mutex));
            this->receiverCondition->notify_one();
//            pthread_cond_signal(&receiverCondition);
            std::unique_lock<std::mutex> lock(*(this->mutex));
            this->workerCondition->wait(lock);
//            pthread_cond_wait(&workerCondition, &(this->mutex));
        }

        if (!isRunning()) {
            //this->mutex->unlock();
            break;
        }
        //this->mutex->unlock();
    }

    printf("Worker %d finished\n", processID);
}

void Worker::executeCurrentTask() {
    while (true) {
        Task task = Task();


        this->mutex->lock();
        if (this->taskQueue->isEmpty()) {
            this->mutex->unlock();
            break;
        }
        task = this->taskQueue->pop();
        this->mutex->unlock();

        printf("Worker %d executing task %d of process %d and weight %d\n",
                processID,
                task.getID(),
                task.getProcessID(),
                task.getWeight());
        usleep(task.getWeight());
    }
}

bool Worker::isRunning() const {
    return this->running;
}

void Worker::initTasks() {
    const int sumWeight = 50000000;
    std::cout << this->taskCount << " " << this->processCount << std::endl;
    int minWeight = 2 * sumWeight / (this->taskCount * (this->processCount + 1));
    int id = 0;
    for (int i = 0; i < this->taskCount; ++i) {
        int weight = minWeight * (i % this->processCount + 1);
        Task task = Task(id, this->processID, weight);
        std::cout << i << "'th task: " << task.to_string() << std::endl;
        if (i % this->processCount == this->processID) {
            std::cout << "pushing..." << std::endl;
            this->mutex->lock();
            this->taskQueue->push(task);
            this->mutex->unlock();
            ++id;
            this->currentProcessSumWeight += task.getWeight();
            std::cout << "pushed." << std::endl;
        }
    }
}
