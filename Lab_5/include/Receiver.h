#ifndef LAB_5_RECEIVER_H
#define LAB_5_RECEIVER_H


#include "TaskQueue.h"
#include "Worker.h"
#include "Sender.h"
#include <condition_variable>
#include <set>

class Receiver {
private:
    const int processID;
    const int processCount;
    std::shared_ptr<TaskQueue> const taskQueue;
    std::shared_ptr<std::mutex> const mutex;
    bool running;

    std::shared_ptr<std::condition_variable> const workerCondition;
    std::shared_ptr<std::condition_variable> const receiverCondition;

    std::shared_ptr<Worker> const worker;
    std::shared_ptr<Sender> const sender;

    static const int taskReplyMPITag = 0xab;
    static const int taskCountRequestMPITag = 0xba;
    static const int taskCountReplyMPITag = 0xbb;
    static const int endingSignal = 404;

    static const int maxReceivingTasksCount = 5;
    std::set<int> otherProcessesWithTasks;

    [[nodiscard]] bool isRunning() const;
    void stop();
public:
    Receiver(int processID,
             int processCount,
             std::shared_ptr<TaskQueue> taskQueue,
             std::shared_ptr<std::mutex> mutex,
             std::shared_ptr<std::condition_variable> workerCondition,
             std::shared_ptr<std::condition_variable> receiverCondition,
             std::shared_ptr<Worker> worker,
             std::shared_ptr<Sender> sender
    );
    void start();
    [[nodiscard]] std::string to_string() const;
};


#endif //LAB_5_RECEIVER_H
