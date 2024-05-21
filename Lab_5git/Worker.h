#ifndef LAB_5_WORKER_H
#define LAB_5_WORKER_H


#include "TaskQueue.h"
#include "Defines.h"
#include <condition_variable>

class Worker {
private:
    const int processID;
    const int processCount;
    TaskQueue* taskQueue;
    std::mutex* mutex;
    bool running;

    std::condition_variable_any* workerCondition;
    std::condition_variable_any* receiverCondition;

    int taskCount;
    int startSumWeight;
    int endSumWeight;
    const int totalSumWeight;// = 47000000;


    void initTasks();
    void executeCurrentTask();

    [[nodiscard]] bool isRunning() const;

public:
    explicit Worker(int processID,
                    int processCount,
                    TaskQueue* taskQueue,
                    std::mutex* mutex,
                    std::condition_variable_any* workerCondition,
                    std::condition_variable_any* receiverCondition,
                    int taskCount,
                    int totalSumWeight
    );
    void start();
    void stop();
    [[nodiscard]] int getStartSumWeight() const;
    [[nodiscard]] int getEndSumWeight() const;
    [[nodiscard]] int getProcessID() const;
    [[nodiscard]] int getProcessCount() const;
};


#endif //LAB_5_WORKER_H
