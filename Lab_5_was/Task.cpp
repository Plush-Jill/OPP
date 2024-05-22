#include "Task.h"

Task::Task(int ID, int processID, int weight) {
    this->ID = ID;
    this->processID = processID;
    this->weight = weight;

}
Task Task::createEmptyTask(int processID) {
    Task task = Task(EMPTY_QUEUE_RESPONSE, processID, 0);
    return task;
}

int Task::getID() const {
    return ID;
}
int Task::getProcessID() const {
    return processID;
}
int Task::getWeight() const {
    return weight;
}

bool Task::isEmpty() const {
    return this->ID == EMPTY_QUEUE_RESPONSE && this->weight == 0;
}

std::string Task::to_string() {
    std::string string = "[TasK: ";
    string += std::to_string(this->ID) + ", ";
    string += std::to_string(this->processID) + ", ";
    string += std::to_string(this->weight) + "]";

    return string;
}

