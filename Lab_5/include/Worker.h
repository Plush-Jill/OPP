#ifndef LAB_5_WORKER_H
#define LAB_5_WORKER_H


#include "TaskQueue.h"
#include <condition_variable>

class Worker {
private:
    const int processID;
    const int processCount;
    std::shared_ptr<TaskQueue> const taskQueue;
    std::shared_ptr<std::mutex> const mutex;
    bool running;

    std::shared_ptr<std::condition_variable> const workerCondition;
    std::shared_ptr<std::condition_variable> const receiverCondition;

    int totalTasksCount;
    int thisWorkerStartCount;
    int startSumWeight;
    int endSumWeight;
    const int totalSumWeight;
    static const int timeScaleToExecuteEachTask = 5'000;
    double sumForAvoidingCompilerOptimization;
    const int taskCountLimitBeforeNotifyingReceiver = 5;

    void initTasks();
    void executeCurrentTask();

    [[nodiscard]] bool isRunning() const;
    [[nodiscard]] bool isEnoughTasksRemainsForDoOnlyOwnTasks() const;

public:
    explicit Worker(int processID,
                    int processCount,
                    std::shared_ptr<TaskQueue> taskQueue,
                    std::shared_ptr<std::mutex> mutex,
                    std::shared_ptr<std::condition_variable> workerCondition,
                    std::shared_ptr<std::condition_variable> receiverCondition,
                    int totalTasksCount,
                    int totalSumWeight
    );
    void start();
    void stop();
    [[nodiscard]] int getStartSumWeight() const;
    [[nodiscard]] int getEndSumWeight() const;
    [[nodiscard]] int getProcessID() const;
    [[nodiscard]] int getProcessCount() const;
    [[nodiscard]] std::string to_string() const;
    [[nodiscard]] double getPartFromTotalSum() const;

    void startWithoutBalancing();
};


#endif //LAB_5_WORKER_H
