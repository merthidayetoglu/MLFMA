#include "vars.h"

//TIMES
double timetotal;
double timetemp;

//GEOMETRY
int numnode;
double *nodepos;
double leafbox;
int numlevel;

//TOPOLOGY
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

  //READ INPUT PARAMETERS
  char *chartemp;
  chartemp = getenv("LEAFBOX");
  double leafbox = atof(chartemp);
  char *modelfile = getenv("MODELFILE");
  chartemp = getenv("SCALE");
  double scalefactor = atof(chartemp);

  //REPORT INPUT PARAMETERS
  if(myrank == 0){
    printf("\n");
    printf("LEAF-LEVEL BOX SIZE: %f\n",leafbox);
    printf("MODEL FILE: %s\n",modelfile);
    printf("SCALE FACTOR: %f\n",scalefactor);
    printf("\n");
    printf("NUMBER OF PROCESSES: %d\n",numrank);
    printf("NUMBER OF THREADS: %d\n",numthread);
    printf("\n");
  }

  //READ GEOMETRY
  MPI_Barrier(MPI_COMM_WORLD);
  timetemp = MPI_Wtime();
  char char_temp[80];
  FILE *pFile = fopen(modelfile,"r");
  while(!feof(pFile)){
    fscanf(pFile,"%s",char_temp);
    if(!strcmp(char_temp,"POINTS")){
      fscanf(pFile,"%s",char_temp);
      numnode = atoi(char_temp);
      fscanf(pFile,"%s",char_temp);
      break;
    }
  }
  if(myrank == 0)printf("NUMBER OF NODES: %d (%f GB)\n",numnode,numnode*sizeof(double)*3/1.0e9);
  nodepos = new double[numnode*3];
  for(int m = 0; m < numnode; m++)
    for(int n = 0; n < 3; n++)
      fscanf(pFile,"%le\n",nodepos+m*3+n);
  fclose(pFile);
  //SCALE GEOMETRY
  #pragma omp parallel for
  for(int n = 0; n < numnode*3; n++)
    nodepos[n] *= scalefactor;
  //#pragma omp parallel for
  //for(int n = 0; n < numnode; n++)
  //  nodepos[n*3] += 5;
  MPI_Barrier(MPI_COMM_WORLD);
  if(myrank == 0)printf("INPUT TIME: %e\n",MPI_Wtime()-timetemp);

  //GEOMETRY SIZE
  double limits[3][2];
  limits[0][0] = nodepos[0];
  limits[0][1] = nodepos[0];
  limits[1][0] = nodepos[1];
  limits[1][1] = nodepos[1];
  limits[2][0] = nodepos[2];
  limits[2][1] = nodepos[2];
  for(int n = 0; n < numnode; n++){
    double x = nodepos[n*3+0];
    double y = nodepos[n*3+1];
    double z = nodepos[n*3+2];
    if(x < limits[0][0]) limits[0][0] = x;
    if(x > limits[0][1]) limits[0][1] = x;
    if(y < limits[1][0]) limits[1][0] = y;
    if(y > limits[1][1]) limits[1][1] = y;
    if(z < limits[2][0]) limits[2][0] = z;
    if(z > limits[2][1]) limits[2][1] = z;
  }
  //PRINT LIMITS
  if(myrank == 0){
    printf("MIN MAX X: %f %f SIZE: %f WAVELENGTHS\n",limits[0][0],limits[0][1],limits[0][1]-limits[0][0]);
    printf("MIN MAX Y: %f %f SIZE: %f WAVELENGTHS\n",limits[1][0],limits[1][1],limits[1][1]-limits[1][0]);
    printf("MIN MAX Z: %f %f SIZE: %f WAVELENGTHS\n",limits[2][0],limits[2][1],limits[2][1]-limits[2][0]);
    printf("GEOMETRIC CENTER: (%f %f %f)\n",(limits[0][1]+limits[0][0])/2,(limits[1][1]+limits[1][0])/2,(limits[2][1]+limits[2][0])/2);
    printf("\n");
  }
  //CONSTRUCT TREE STRUCTURE
  double max = 0;
  for(int m = 0; m < 3; m++)
    for(int n = 0; n < 2; n++)
      if(limits[m][n] > max)
        max = abs(limits[m][n]);
  numlevel = 2;
  double size = leafbox;
  while(size < max){
    numlevel++;
    size *= 2;
  }
  double boxsize[numlevel];
  boxsize[numlevel-1] = leafbox;
  for(int l = numlevel-1; l > 0; l--)
    boxsize[l-1] = boxsize[l]*2;
  if(myrank == 0){
    printf("NUMBER OF LEVELS: %d\n",numlevel);
    printf("LARGEST BOX SIZE: %f\n",size*2);
    printf("FULL TREE STRUCTURE\n");
    long allbox = 1;
    long box = 8;
    for(int l = 1; l < numlevel; l++){
      printf("level: %d numbox: %ld boxsize: %f\n",l,box,boxsize[l]);
      allbox += box;
      box *= 8;
    }
    printf("TOTAL NUMBER OF BOXES: %ld\n\n",allbox);
    printf("PRUNING TREE STRUCTURE\n");
  }
  long numbox[numlevel];
  numbox[0] = 1;
  double *boxcenter[numlevel];
  boxcenter[0] = new double[numbox[0]*3];
  boxcenter[0][0] = 0.0;
  boxcenter[0][1] = 0.0;
  boxcenter[0][2] = 0.0;
  long *nodedispl[numlevel];
  nodedispl[0] = new long[numbox[0]+1];
  nodedispl[0][0] = 0;
  nodedispl[0][1] = numnode;
  long *childdispl[numlevel-1];
  childdispl[0] = new long[numbox[0]+1];
  childdispl[0][0] = 0;

  char *childtag = new char[numnode];
  long *sortindex = new long[numnode];
  for(long n = 0; n < numnode; n++)
    sortindex[n] = n;

  if(myrank == 0)
    for(int l = 0; l < 1; l++){
      printf("LEVEL %d\n",l);
      //FIND NUMBER OF CHILDREN FOR EACH PARENT BOX
      for(int parent = 0; parent < numbox[l]; parent++){
        double cenx = boxcenter[l][parent*3+0];
        double ceny = boxcenter[l][parent*3+1];
        double cenz = boxcenter[l][parent*3+2];
        double corx = cenx-boxsize[l]/2;
        double cory = ceny-boxsize[l]/2;
        double corz = cenz-boxsize[l]/2;
        printf("parent: %d center(%f %f %f) corner(%f %f %f)\n",parent,cenx,ceny,cenz,corx,cory,corz);
	long bin[8] = {0};
	for(long n = nodedispl[l][parent]; n < nodedispl[l][parent+1]; n++){
          long child = sortindex[n];
          int binx = (int)((nodepos[child*3+0]-corx)/boxsize[l+1]);
          int biny = (int)((nodepos[child*3+1]-cory)/boxsize[l+1]);
          int binz = (int)((nodepos[child*3+2]-corz)/boxsize[l+1]);
	  int binind = 4*binz+2*biny+binx;
	  bin[binind]++;
	  childtag[n] = binind;
	}
	printf("bin: ");
	int numchild = 0;
	for(int b = 0; b < 8; b++){
          if(bin[b]){
            printf("%ld ",bin[b]);
	    bin[b] = numchild;
	    numchild++;
	  }
	}
        printf("\n");
        printf("number of children: %d\n",numchild);
        childdispl[l][parent+1] = childdispl[l][parent] + numchild;
      }
      printf("\n");
    }

  MPI_Barrier(MPI_COMM_WORLD);
  if(myrank == 0)printf("\nTOTAL TIME: %e\n\n",MPI_Wtime()-timetotal);

  MPI_Finalize();

  return 0;
}
