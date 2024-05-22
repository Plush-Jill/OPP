#include "../include/TaskQueue.h"




bool TaskQueue::isEmpty() const {
    return this->queue->empty();
}
void TaskQueue::push(Task task) {
    this->queue->push(task);
}
Task TaskQueue::pop() {
    Task task = this->queue->front();
    this->queue->pop();
    return task;
}
TaskQueue::TaskQueue(int capacity) :
queue(std::make_shared<std::queue<Task>>())
{
}
int TaskQueue::getSize() const {
    return static_cast<int>(this->queue->size());
}