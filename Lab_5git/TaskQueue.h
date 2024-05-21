#ifndef LAB_5_TASKQUEUE_H
#define LAB_5_TASKQUEUE_H


#include <memory>
#include <vector>
#include <queue>
#include "Task.h"

class TaskQueue {
private:
    std::queue<Task>* queue;
    int capacity;

public:
    explicit TaskQueue(int capacity);
    explicit TaskQueue();

    [[nodiscard]] bool isEmpty() const;
    [[nodiscard]] bool isFull() const;
    void push(Task task);
    Task pop();

    [[nodiscard]] int getCapacity() const;
    [[nodiscard]] int getSize() const;

};

#endif //LAB_5_TASKQUEUE_H
