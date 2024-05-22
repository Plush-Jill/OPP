#include <iostream>
#include "TaskQueue.h"



TaskQueue::TaskQueue(int capacity) : TaskQueue(){
    this->capacity = capacity;
}
TaskQueue::TaskQueue() {
    this->tasks = std::make_shared<std::vector<Task>>();
    this->tasks->reserve(capacity);
    this->capacity = 2000;
    this->count = 0;
    this->popIndex = 0;
}

bool TaskQueue::isEmpty() const {
    return this->count == 0;
}

bool TaskQueue::isFull() const {
    return this->count == this->capacity;
}

void TaskQueue::push(Task task) {
    std::cout << "pushing " << 1 << std::endl;
    /*if (isFull()){
        throw std::exception();
    }*/
    std::cout << "pushing " << 2 << std::endl;

    int pushIndex = (this->popIndex + this->count) % this->capacity;
    std::cout << "pushing " << 3 << std::endl;

    this->tasks->assign(pushIndex, task);// = task;
    std::cout << "pushing " << 4 << std::endl;

    ++this->count;
    std::cout << "pushing " << 5 << std::endl;

}
const std::shared_ptr<std::vector<Task>> &TaskQueue::getTasks() const {
    return tasks;
}
int TaskQueue::getCapacity() const {
    return capacity;
}
int TaskQueue::getCount() const {
    return count;
}
int TaskQueue::getPopIndex() const {
    return popIndex;
}

Task TaskQueue::pop() {
    Task task = this->tasks->at(this->popIndex);
    this->popIndex = (this->popIndex + 1) & this->capacity;
    --this->count;

    return task;
}
