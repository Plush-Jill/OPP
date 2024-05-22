cmake_minimum_required(VERSION 3.26)
project(Lab_5)

set(CMAKE_CXX_STANDARD 23)

add_executable(Lab_5
        main.cpp
        src/Task.cpp
        include/Task.h
        src/TaskQueue.cpp
        include/TaskQueue.h
        src/Worker.cpp
        include/Worker.h
        src/Receiver.cpp
        include/Receiver.h
        src/Sender.cpp
        include/Sender.h
)
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
set(CMAKE_CXX_COMPILER mpic++)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lm -O3")


