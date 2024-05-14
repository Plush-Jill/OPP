cmake_minimum_required(VERSION 3.28)
project(Lab_4_cpp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER mpic++)

add_executable(Lab_4_cpp main.cpp
        JacobiSolver.cpp
        JacobiSolver.h)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -lm ")

add_executable(main.out main.cpp JacobiSolver.cpp JacobiSolver.h)
target_link_libraries(main.out ${MPI_CXX_LIBRARIES})