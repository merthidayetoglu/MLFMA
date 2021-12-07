#!/bin/bash

date

export LEAFBOX=0.25
export MODELFILE=sphere_0p01m.vtk
export SCALE=20

#mpirun -n 1 -N 1 -x OMP_NUM_THREADS=40 -x OMP_PLACES=cores ./MLFMA
#jsrun -n1 -a1 -g1 -c21 -EOMP_NUM_THREADS=21 -r1 -bpacked:21 js_task_info ./MLFMA
mpirun -n 2 -N 2 --map-by node:PE=8 -x OMP_NUM_THREADS=8 -x OMP_PLACES=cores ./MLFMA

date
