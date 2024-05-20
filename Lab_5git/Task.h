#ifndef LAB_5_TASK_H
#define LAB_5_TASK_H


#include <string>

class Task {
private:
    int ID;
    int processID;
    int weight;

    static const int EMPTY_QUEUE_RESPONSE = -1;
public:
    Task() = default;
    Task(int ID, int processID, int weight);
    static Task createEmptyTask(int processID);

    [[nodiscard]] int getID() const;
    [[nodiscard]] int getProcessID() const;
    [[nodiscard]] int getWeight() const;
    [[nodiscard]] bool isEmpty() const;

    std::string to_string();
};


#endif //LAB_5_TASK_H
