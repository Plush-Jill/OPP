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

    std::condition_variable* workerCondition;
    std::condition_variable* receiverCondition;

    int taskCount;
    int currentProcessStartSumWeight;
    int currentProcessEndSumWeight;
    const int totalSumWeight;// = 47000000;

    pthread_mutex_t* mutexC;
    pthread_cond_t* workerConditionC;
    pthread_cond_t* receiverConditionC;

    void initTasks();
    void executeCurrentTask();

    [[nodiscard]] bool isRunning() const;

public:
    explicit Worker(int processID,
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
    );
    void start();
    void stop();
};


#endif //LAB_5_WORKER_H
