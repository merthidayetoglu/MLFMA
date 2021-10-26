#include <stdio.h>
#include "mpi.h"
#include "omp.h"

int myrank;
int numrank;

int main(int argc, char** argv){

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  MPI_Comm_size(MPI_COMM_WORLD,&numrank);

  MPI_Finalize();

  return 0;
}
