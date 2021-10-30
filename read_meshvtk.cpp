#include <fstream>
#include <stdio.h>
#include <string.h>

using namespace std;

extern char modelname[];
extern int num_nodes;
extern int num_tris;
extern double *node_pos;
extern int *tri_nodes;

void read_meshvtk(){
    //
    //
    // read_meshvtk
    // Reading the node coordinates and triangle nodes from vtk file.
    //
    //
    char char_temp[80];
    string string_temp;
    FILE *pFile;

    pFile = fopen (modelname,"r");
    if (pFile==NULL){
        printf("Model File Not Found! \n");
        exit(EXIT_FAILURE);
    } else{
        printf("Model Node File: %s \n", modelname);
    }

    //
    //Reading node positions
    //
    while(!feof(pFile)){
        fscanf(pFile, "%s", char_temp);
        if(!strcmp(char_temp,"POINTS")){
            fscanf(pFile, "%s", char_temp);
            num_nodes=stoi(char_temp);
            printf("Number of Nodes: %d\n", num_nodes);
            fscanf(pFile, "%s", char_temp);

            break;
        }
    }
    node_pos = new double[num_nodes*3];
    for(int indline=0;indline<num_nodes;indline++){
        fscanf (pFile, "%lf %lf %lf\n", &node_pos[indline*3+0], &node_pos[indline*3+1], &node_pos[indline*3+2]);
    }

    //
    // Reading triangle nodes
    //
    while(!feof(pFile)){
        fscanf(pFile, "%s", char_temp);
        if(!strcmp(char_temp,"CELLS")){
            fscanf(pFile, "%s", char_temp);
            num_tris=stoi(char_temp);
            printf("Number of Triangles: %d\n", num_tris);
            fscanf(pFile, "%s", char_temp);
            break;
        }
    }

    tri_nodes = new int[num_tris*3];
    for(int indline=0;indline<num_tris;indline++){
        fscanf(pFile, "%s", char_temp);
        fscanf(pFile, "%d %d %d\n", &tri_nodes[indline*3+0], &tri_nodes[indline*3+1], &tri_nodes[indline*3+2]);
    }

    fclose (pFile);

}
