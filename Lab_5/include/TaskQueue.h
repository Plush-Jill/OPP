#ifndef LAB_5_TASKQUEUE_H
#define LAB_5_TASKQUEUE_H


#include <memory>
#include <vector>
#include <queue>
#include "Task.h"

class TaskQueue {
private:
    std::shared_ptr<std::queue<Task>> const queue;

public:
    explicit TaskQueue(int capacity);

    [[nodiscard]] bool isEmpty() const;
    void push(Task task);
    Task pop();

    [[nodiscard]] int getRemainsTasksCount() const;
};

#endif //LAB_5_TASKQUEUE_H
