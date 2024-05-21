#include <iostream>
#include "TaskQueue.h"




bool TaskQueue::isEmpty() const {
    return this->queue->empty();
}

bool TaskQueue::isFull() const {
    return getSize() == this->capacity;
}

void TaskQueue::push(Task task) {
    this->queue->push(task);
}
int TaskQueue::getCapacity() const {
    return capacity;
}
Task TaskQueue::pop() {
    Task task = this->queue->front();
    this->queue->pop();
    return task;
    //return this->queue->front();
}

TaskQueue::TaskQueue(int capacity) {
    this->queue = new std::queue<Task>();
    this->capacity = capacity;
}

int TaskQueue::getSize() const {
    return static_cast<int>(this->queue->size());
}

TaskQueue::TaskQueue() : TaskQueue(4000){

}
