#include "../include/Task.h"

Task::Task(int ID,
           int processID,
           int weight) :
           ID(ID), processID(processID), weight(weight)
           {

           }
Task Task::createEmptyTask(int processID) {
    Task task = Task(0, processID, 0);
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
    return this->ID == 0 && this->weight == 0;
}
std::string Task::to_string() const {
    std::string string = "[Task: ";
    string += "ID: " + std::to_string(this->ID) + ", ";
    string += "PID: " + std::to_string(this->processID) + ", ";
    string += "W: " + std::to_string(this->weight) + "]";

    return string;
}

Task::Task() : ID(0), processID(0), weight(0){

}

