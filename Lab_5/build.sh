#!/bin/bash



module load mpi/mpich-x86_64

/home/fit_opp/22202/mpe2-2.4.9b/bin/mpecxx main.cpp -c -pthread -O3 -lm
/home/fit_opp/22202/mpe2-2.4.9b/bin/mpecxx Task.cpp -c -pthread -O3 -lm
/home/fit_opp/22202/mpe2-2.4.9b/bin/mpecxx TaskQueue.cpp -c -pthread -O3 -lm
/home/fit_opp/22202/mpe2-2.4.9b/bin/mpecxx Worker.cpp -c -pthread -O3 -lm
/home/fit_opp/22202/mpe2-2.4.9b/bin/mpecxx Receiver.cpp -c -pthread -O3 -lm
/home/fit_opp/22202/mpe2-2.4.9b/bin/mpecxx Sender.cpp -c -pthread -O3 -lm


/home/fit_opp/22202/mpe2-2.4.9b/bin/mpecxx main.o Task.o TaskQueue.o Worker.o Receiver.o Sender.o -pthread -o main








module unload mpi/mpich-x86_64
