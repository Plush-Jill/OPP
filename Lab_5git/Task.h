#ifndef LAB_5_TASK_H
#define LAB_5_TASK_H


#include <string>

class Task {
private:
    int ID;
    int processID;
    int weight;

public:
    Task();
    Task(int ID, int processID, int weight);
    static Task createEmptyTask(int processID);

    [[nodiscard]] int getID() const;
    [[nodiscard]] int getProcessID() const;
    [[nodiscard]] int getWeight() const;
    [[nodiscard]] bool isEmpty() const;

    [[nodiscard]] std::string to_string() const;
};


#endif //LAB_5_TASK_H
