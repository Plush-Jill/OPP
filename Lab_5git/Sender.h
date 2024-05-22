#ifndef LAB_5_SENDER_H
#define LAB_5_SENDER_H


#include <condition_variable>
#include "TaskQueue.h"

class Sender {
private:
    const int processID;
    const int processCount;
    std::shared_ptr<TaskQueue> const taskQueue;
    std::shared_ptr<std::mutex> const mutex;
    bool running;

    std::shared_ptr<std::condition_variable> const workerCondition;
    std::shared_ptr<std::condition_variable> const receiverCondition;

    static const int taskRequestMPITag = 0xa;
    static const int taskReplyMPITag = 0xb;
    static const int endingSignal = 404;

public:
    explicit Sender(int processID,
                    int processCount,
                    std::shared_ptr<TaskQueue> taskQueue,
                    std::shared_ptr<std::mutex> mutex,
                    std::shared_ptr<std::condition_variable> workerCondition,
                    std::shared_ptr<std::condition_variable> receiverCondition
    );
    void start();
    void stop();
};


#endif //LAB_5_SENDER_H
