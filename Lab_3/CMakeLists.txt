cmake_minimum_required(VERSION 3.28)
project(Lab_3 C)

set(CMAKE_C_STANDARD 23)
find_package(MPI REQUIRED)
add_executable(Lab_3
        main.c)
target_include_directories(Lab_3
        PUBLIC ${MPI_C_INCLUDE_DIRS})

