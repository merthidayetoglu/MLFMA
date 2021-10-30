#include <stdio.h>
#include <complex>

using namespace std;

// Variables
// created in read_simulationinput
char modelname[200];
int num_mlfmalevels;
double val_freq;

// created in read_meshmodel
int num_nodes;
int num_tris;
double *node_pos;
int *tri_nodes;

// created in genetare_constants
double pi=3.1415926535897932384626433;
double val_e=2.7182818284590452353602874;
double val_epsilonO=8.854187817620390e-12;
double val_muO=1.2566370621219e-6;
double val_omega;
double val_c;
double val_lambdaO;
double val_kO;
double val_etaO;

// created in generate_octtree
int *num_truncation;
int *num_truncation_pc;
double *size_boxes;
int num_filledboxesall;
int *num_filledboxes;
int *num_children;
double *box_pos;
int *children;
int *children2parent;
int *rwg_ind;


// Functions
void read_simulationinput();
void read_meshmodel();
void genetare_constants();
void generate_octtree();

