#ifndef LAB_5_TASKQUEUE_H
#define LAB_5_TASKQUEUE_H


#include <memory>
#include <vector>
#include <queue>
#include "Task.h"

class TaskQueue {
private:
    std::shared_ptr<std::vector<Task>> tasks;
    std::queue<Task> queue;
    int capacity;
    int count;
    int popIndex;

public:
    explicit TaskQueue(int capacity);
    explicit TaskQueue();

    [[nodiscard]] bool isEmpty() const;
    [[nodiscard]] bool isFull() const;
    void push(Task task);
    Task pop();

    [[nodiscard]] const std::shared_ptr<std::vector<Task>> &getTasks() const;
    [[nodiscard]] int getCapacity() const;
    [[nodiscard]] int getCount() const;
    [[nodiscard]] int getPopIndex() const;

};

#endif //LAB_5_TASKQUEUE_H
