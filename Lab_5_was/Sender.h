#ifndef LAB_5_SENDER_H
#define LAB_5_SENDER_H


#include "TaskQueue.h"
#include "Defines.h"

class Sender {
private:
    int processID;
    std::shared_ptr<TaskQueue> taskQueue;
    std::shared_ptr<std::mutex> mutex;
    bool running;

public:
    Sender(int processID, std::shared_ptr<TaskQueue>& taskQueue, std::shared_ptr<std::mutex>& mutex);
    void start();
};


#endif //LAB_5_SENDER_H
