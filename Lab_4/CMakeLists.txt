cmake_minimum_required(VERSION 3.28)
project(Lab_4)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER mpic++)

add_executable(Lab_4 main.cpp
        JacobiSolver.cpp
        JacobiSolver.h)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -lm ")

