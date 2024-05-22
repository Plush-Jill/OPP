#!/bin/bash

mpic++ -c -std=c++20 -O3 -lm main.cpp
mpic++ -c -std=c++20 -O3 -lm Task.cpp
#mpic++ -c -std=c++20 -O3 -lm Task.h
mpic++ -c -std=c++20 -O3 -lm TaskQueue.cpp
#mpic++ -c -std=c++20 -O3 -lm TaskQueue.h
#mpic++ -c -std=c++20 -O3 -lm Exceptions.h
mpic++ -c -std=c++20 -O3 -lm Worker.cpp
#mpic++ -c -std=c++20 -O3 -lm Worker.h
mpic++ -c -std=c++20 -O3 -lm Receiver.cpp
#mpic++ -c -std=c++20 -O3 -lm Receiver.h
mpic++ -c -std=c++20 -O3 -lm Sender.cpp
#mpic++ -c -std=c++20 -O3 -lm Sender.h
#mpic++ -c -std=c++20 -O3 -lm Defines.h

mpic++ -std=c++20 -fmodules-ts -o main main.o Task.o TaskQueue.o Worker.o Receiver.o Sender.o
# Task_impl.o TaskQueue_impl.o