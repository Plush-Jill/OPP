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
)
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
set(CMAKE_CXX_COMPILER mpic++)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lm -O3")
add_executable(Lab_5c main.c)

#cmake_minimum_required(VERSION 3.20)
#project(Lab_5)
#
#set(CMAKE_CXX_STANDARD 20)
#
#add_executable(Lab_5
#        src/main.cpp
#        src/Task.cpp
#        src/Task.h
#        src/TaskQueue.cpp
#        src/TaskQueue.h
#        src/Exceptions.h
#        src/Worker.cpp
#        src/Worker.h
#        src/Receiver.cpp
#        src/Receiver.h
#        src/Sender.cpp
#        src/Sender.h
#)
##find_package(MPI REQUIRED)
##include_directories(${MPI_INCLUDE_PATH})
#set(CMAKE_CXX_COMPILER /home/fit_opp/22202/mpe2-2.4.9b/bin/mpecxx)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -pthread -lm -O3")

