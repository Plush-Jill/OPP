#include <exception>
#include <string>

class TaskQueueException : std::exception{
private:
    std::string message_;
public:
    explicit TaskQueueException(std::string  message) : message_(std::move(message)) {}
    explicit TaskQueueException() = default;
    virtual const char* what() const noexcept{
        return message_.c_str();
    }
};

class InitMPIException : std::exception{
private:
    std::string message_;
public:
    explicit InitMPIException(std::string  message) : message_(std::move(message)) {}
    explicit InitMPIException() = default;
    virtual const char* what() const noexcept{
        return message_.c_str();
    }
};