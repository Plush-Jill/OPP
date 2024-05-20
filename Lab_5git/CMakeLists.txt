cmake_minimum_required(VERSION 3.26)
project(Lab_5)

set(CMAKE_CXX_STANDARD 23)

add_executable(Lab_5
        main.cpp
        Task.cpp
        Task.h
        TaskQueue.cpp
        TaskQueue.h
        Exceptions.h
        Worker.cpp
        Worker.h
        Receiver.cpp
        Receiver.h
        Sender.cpp
        Sender.h
        Defines.h
)
add_executable(Lab_5c mainC.c)
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
set(CMAKE_CXX_COMPILER mpic++)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lm")
add_executable(n5 n5.cpp)
