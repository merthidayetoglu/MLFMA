#!/bin/bash

date
module load cuda

#mpirun -n 1 -N 1 -x OMP_NUM_THREADS=40 -x OMP_PLACES=cores ./MLFMA
jsrun -n1 -a1 -g1 -c21 -EOMP_NUM_THREADS=21 -r1 -bpacked:21 js_task_info ./MLFMA

date
