#ifndef LAB_5_RECEIVER_H
#define LAB_5_RECEIVER_H


#include "TaskQueue.h"
#include "Defines.h"
#include <condition_variable>

class Receiver {
private:
    int processID;
    std::shared_ptr<TaskQueue> taskQueue;
    std::shared_ptr<std::mutex> mutex;
    bool running;

    int processCount;
    std::shared_ptr<std::condition_variable> workerCondition;
    std::shared_ptr<std::condition_variable> receiverCondition;

    bool isRunning() const;
public:
    Receiver(int processID, int processCount, std::shared_ptr<TaskQueue> &taskQueue, std::shared_ptr<std::mutex> &mutex,
             std::shared_ptr<std::condition_variable> &receiverCondition,
             std::shared_ptr<std::condition_variable> &workerCondition);
    void start();
};


#endif //LAB_5_RECEIVER_H
