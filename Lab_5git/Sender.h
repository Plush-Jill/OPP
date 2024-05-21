#ifndef LAB_5_SENDER_H
#define LAB_5_SENDER_H


#include <condition_variable>
#include "TaskQueue.h"
#include "Defines.h"

class Sender {
private:
    const int processID;
    const int processCount;
    TaskQueue* taskQueue;
    std::mutex* mutex;
    bool running;

    std::condition_variable_any* workerCondition;
    std::condition_variable_any* receiverCondition;

public:
    explicit Sender(int processID,
                    int processCount,
                    TaskQueue* taskQueue,
                    std::mutex* mutex,
                    std::condition_variable_any* workerCondition,
                    std::condition_variable_any* receiverCondition
    );
    void start();
    void stop();
};


#endif //LAB_5_SENDER_H
