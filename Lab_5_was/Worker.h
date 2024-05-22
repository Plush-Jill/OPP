#ifndef LAB_5_WORKER_H
#define LAB_5_WORKER_H


#include "TaskQueue.h"
#include "Defines.h"
#include <condition_variable>

class Worker {
private:
    int processID;
    std::shared_ptr<TaskQueue> taskQueue;
    std::shared_ptr<std::mutex> mutex;
    bool running;

    int taskCount;
    int currentProcessSumWeight;
    int processCount;

    std::shared_ptr<std::condition_variable> workerCondition;
    std::shared_ptr<std::condition_variable> receiverCondition;

    void initTasks();
    void executeCurrentTask();
    void waitForNewTasks();

    bool isRunning() const;

public:
    Worker(int processID, int processCount, int taskCount, std::shared_ptr<std::mutex> &mutex,
           std::shared_ptr<std::condition_variable> &workerCondition,
           std::shared_ptr<std::condition_variable> &receiverCondition, std::shared_ptr<TaskQueue> &taskQueue);
    void start();
};


#endif //LAB_5_WORKER_H
