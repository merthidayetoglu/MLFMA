#!/bin/bash

date

export LEAFBOX=0.25
export MODELFILE=/gpfs/alpine/csc362/scratch/merth/MLFMA/sphere_0p01m.vtk
export SCALE=1.0;

#mpirun -n 1 -N 1 -x OMP_NUM_THREADS=40 -x OMP_PLACES=cores ./MLFMA
jsrun -n1 -a1 -g1 -c21 -EOMP_NUM_THREADS=21 -r1 -bpacked:21 js_task_info ./MLFMA

date
