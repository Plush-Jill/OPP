#ifndef LAB_5_RECEIVER_H
#define LAB_5_RECEIVER_H


#include "TaskQueue.h"
#include "Defines.h"
#include "Worker.h"
#include "Sender.h"
#include <condition_variable>

class Receiver {
private:
    const int processID;
    const int processCount;
    TaskQueue* taskQueue;
    std::mutex* mutex;
    bool running;

    std::condition_variable* workerCondition;
    std::condition_variable* receiverCondition;

    Worker* worker;
    Sender* sender;

    pthread_mutex_t* mutexC;
    pthread_cond_t* workerConditionC;
    pthread_cond_t* receiverConditionC;

    [[nodiscard]] bool isRunning() const;
    void stop();
public:
    Receiver(int processID,
             int processCount,
             TaskQueue* taskQueue,
             std::mutex* mutex,
             std::condition_variable* workerCondition,
             std::condition_variable* receiverCondition,
             Worker* worker,
             Sender* sender,
             pthread_mutex_t* mutexC,
             pthread_cond_t* workerConditionC,
             pthread_cond_t* receiverConditionC
    );
    void start();
};


#endif //LAB_5_RECEIVER_H
