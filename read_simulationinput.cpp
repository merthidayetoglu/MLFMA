#include <stdio.h>
#include <stdlib.h>

using namespace std;

extern char modelname[];
extern int num_mlfmalevels;
extern double val_freq;

void read_simulationinput(){
    //
    //
    // read_simulationinput
    // Reading the simulation parameters from "input.dat"

    char temp1[200];
    FILE *pFile;
    pFile = fopen ("input.dat","r");
    if (pFile==NULL){
        printf("Input File Not Found! \n");
        exit(EXIT_FAILURE);
    }

    fscanf (pFile, "%11c %s\n", temp1, modelname);
    fscanf (pFile, "%10c %lf\n", temp1, &val_freq);
    fscanf (pFile, "%17c %d\n", temp1, &num_mlfmalevels);
    fclose(pFile);


}



