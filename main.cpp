#include <stdio.h>
#include <mpi.h>
#include <omp.h>

double timetotal;
double timetemp;

double leafbox;

int numnode;
double *nodepos;

int myrank;
int numrank;
int numthread;

int main(int argc, char** argv) {

  timetotal = MPI_Wtime();

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  MPI_Comm_size(MPI_COMM_WORLD,&numrank);

  #pragma omp parallel
  if(omp_get_thread_num() == 0)
    numthread = omp_get_num_threads();

  //INPUT PARAMETERS
  char *chartemp;
  chartemp = getenv("LEAFBOX");
  double leafbox = atof(chartemp);
  char *modelfile = getenv("MODELFILE");
  chartemp = getenv("SCALE");
  double scalefactor = atof(chartemp);

  //REPORT INPUT PARAMETERS
  if(myrank == 0){
    printf("LEAF-LEVEL BOX SIZE: %e\n",leafbox);
    printf("MODEL FILE: %s\n",modelfile);
    printf("SCALE FACTOR: %e\n",scalefactor);
    printf("\n");
    printf("NUMBER OF PROCESSES: %d\n",numrank);
    printf("NUMBER OF THREADS: %d\n",numthread);
    printf("\n");
  }

  //READ GEOMETRY
  MPI_Barrier(MPI_COMM_WORLD);
  timetemp = MPI_Wtime();
  if(myrank == 0)printf("MODEL FILE: %s\n",modelfile);
  if(myrank == 0)printf("SCALE: %e\n",scalefactor);
  FILE *pFile = fopen(modelfile,"rb");
  fread(&numnode,sizeof(int),1,pFile);
  if(myrank == 0)printf("NUMBER OF NODES: %d (%f GB)\n",numnode,numnode*sizeof(double)*3/1.0e9);
  nodepos = new double[numnode*3];
  fread(nodepos,sizeof(double),numnode*3,pFile);
  fclose(pFile);
  for(int n = 0; n < numnode*3; n++)
    nodepos[n] *= scalefactor;
  MPI_Barrier(MPI_COMM_WORLD);
  if(myrank == 0)printf("INPUT TIME: %e\n",MPI_Wtime()-timetemp);

  //CONSTRUCT TREE STRUCTURE

  MPI_Barrier(MPI_COMM_WORLD);
  if(myrank == 0)printf("TOTAL TIME: %e\n",MPI_Wtime()-timetotal);

  MPI_Finalize();

  return 0;
}
