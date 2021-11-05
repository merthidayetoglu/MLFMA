#include <stdio.h>
#include <stdlib.h>
#include <cmath>

using namespace std;

extern int num_nodes;
extern int num_tris;
extern double *node_pos;
extern int *tri_nodes;
extern int *rwg_nodes;
extern double val_lambdaO;
extern double val_kO;

extern int num_mlfmalevels;
extern double *size_boxes;
extern int *num_truncation;

extern int num_digits;

extern int num_filledboxesall;
extern int *num_filledboxes;
extern int *num_children;
extern double *box_pos;
extern int *children;
extern int *children2parent;

struct sortedElements {
    int elemID;
    long int elemVal;
};

int compare (const void *a, const void *b);
int compareind (const void *a, const void *b);
long int dec2bin2dec(int,int,int,int);


void generate_octtree(){
    //
    //
    // generate_octtree
    // Construct the oct-tree structure of MLFMA. Determine the clusters/boxes of RWG functions.
    double objectmin[3];
    double objectmax[3];
    double objectcenter[3];
    double objectsize;
    double boxmin[3];
    int nodeposTemp[3];
    int *node_box;
    int temp_nodes[2],temp_nodes_c[3];
    sortedElements *node_box_forsort;


    // Determine the min-max location, the center point, and the object size of the mesh-geometry.
    objectmin[0]=0.0;
    objectmin[1]=0.0;
    objectmin[2]=0.0;
    objectmax[0]=0.0;
    objectmax[1]=0.0;
    objectmax[2]=0.0;
    for(int ind1=0; ind1<num_nodes; ind1++){
        if(objectmin[0]>node_pos[ind1*3+0])
            objectmin[0]=node_pos[ind1*3+0];
        if(objectmin[1]>node_pos[ind1*3+1])
            objectmin[1]=node_pos[ind1*3+1];
        if(objectmin[2]>node_pos[ind1*3+2])
            objectmin[2]=node_pos[ind1*3+2];
        if(objectmax[0]<node_pos[ind1*3+0])
            objectmax[0]=node_pos[ind1*3+0];
        if(objectmax[1]<node_pos[ind1*3+1])
            objectmax[1]=node_pos[ind1*3+1];
        if(objectmax[2]<node_pos[ind1*3+2])
            objectmax[2]=node_pos[ind1*3+2];
    }
    objectcenter[0]=(objectmax[0]+objectmin[0])/2;
    objectcenter[1]=(objectmax[1]+objectmin[1])/2;
    objectcenter[2]=(objectmax[2]+objectmin[2])/2;
    objectsize=0.0;
    for(int ind1=0; ind1<3; ind1++){
        if(objectsize<(objectmax[ind1]-objectmin[ind1]))
            objectsize=objectmax[ind1]-objectmin[ind1];
    }
    printf("Min: %f %f %f \n", objectmin[0],objectmin[1],objectmin[2]);
    printf("Max: %f %f %f \n", objectmax[0],objectmax[1],objectmax[2]);
    printf("Center: %f %f %f \n", objectcenter[0],objectcenter[1],objectcenter[2]);
    printf("Model Size: %f \n", objectsize);
    printf("\n");

    // Determine the number of MLFMA levels for the simulation.
    if(num_mlfmalevels==1){
        num_mlfmalevels=floor(log2(objectsize/val_lambdaO))+3;
        if(num_mlfmalevels<2){
                num_mlfmalevels=2;
        }
    }
    printf("MLFMA Levels: %d \n", num_mlfmalevels);

    // Calculate the required box sizes at each MLFMA level.
    size_boxes=new double[num_mlfmalevels];
    size_boxes[0]=objectsize;
    for(int ind1=1; ind1<num_mlfmalevels; ind1++)
        size_boxes[ind1]=size_boxes[ind1-1]/2;
    double smallestboxsize=size_boxes[num_mlfmalevels-1];
    boxmin[0]=objectcenter[0]-objectsize/2;
    boxmin[1]=objectcenter[1]-objectsize/2;
    boxmin[2]=objectcenter[2]-objectsize/2;

    // Calculate the required truncation numbers (number of harmonics) at each MLFMA level.
    num_truncation=new int[num_mlfmalevels];
    for(int ind1=0; ind1<num_mlfmalevels; ind1++){
        num_truncation[ind1]=ceil((1.73*abs(val_kO)*size_boxes[ind1]+2.16*(pow(num_digits,0.666666666666666))*pow(abs(val_kO)*size_boxes[ind1],0.333333333333333)));
        printf("Level: %d, Num. Harmonics: %d \n", ind1, num_truncation[ind1]);
    }

    // Calculate Node location indices at the least level boxes. Storing the variables in XYZ coordinates and in combined-form.
    // The Node box location indices are sorted for the calculation of box locations at the higher levels.
    node_box_forsort=new sortedElements[num_nodes];
    node_box=new int[num_nodes*3];

    for(int ind1=0; ind1<num_nodes; ind1++){
        node_box_forsort[ind1].elemID=ind1;
        node_box_forsort[ind1].elemVal=0;
    }
    for(int ind1=0; ind1<num_nodes; ind1++){
        node_box[ind1*3+0]=floor((node_pos[ind1*3+0]-boxmin[0])/smallestboxsize);
        node_box[ind1*3+1]=floor((node_pos[ind1*3+0]-boxmin[1])/smallestboxsize);
        node_box[ind1*3+2]=floor((node_pos[ind1*3+0]-boxmin[2])/smallestboxsize);
        node_box_forsort[ind1].elemVal=dec2bin2dec(node_box[ind1*3+0],node_box[ind1*3+1],node_box[ind1*3+2],num_mlfmalevels);
    }
    qsort(node_box_forsort, num_nodes, sizeof(struct sortedElements), compare);


    // Determine the total number of filled boxes "ind_filledboxes"
    // for the allocation of arrays "num_children" and "box_pos"
    num_filledboxes=new int[num_mlfmalevels*3];
	num_filledboxesall=0;
    for(int indlev=num_mlfmalevels; indlev>0; indlev--){
        num_filledboxes[indlev]=1;
        temp_nodes[0]=node_box_forsort[0].elemVal/pow(8,num_mlfmalevels-indlev);
        temp_nodes[1]=node_box_forsort[1].elemVal/pow(8,num_mlfmalevels-indlev);
        for(int indnode=1; indnode<num_nodes; indnode++){
            if (temp_nodes[1]!=temp_nodes[0]){
                num_filledboxes[indlev]++;
            }
            temp_nodes[indnode-1]=node_box_forsort[indnode-1].elemVal/pow(8,num_mlfmalevels-indlev);
            temp_nodes[indnode]=node_box_forsort[indnode].elemVal/pow(8,num_mlfmalevels-indlev);
        }
        num_filledboxesall+=num_filledboxes[indlev];
    }

    // Calculate the box positions in cartesian coordinates "box_pos" and determine the number of children at each box "num_children"
    num_children=new int[num_filledboxesall];
    box_pos=new double[num_filledboxesall*3];
    double boxsizesvec=0.0;
    int ind_filledboxes=-1;
    for(int indlev=num_mlfmalevels-1; indlev>=0; indlev--){
        boxsizesvec=size_boxes[indlev]/2;
        num_filledboxes[indlev]=1;
        num_children[ind_filledboxes+num_filledboxes[indlev]]=1;
        int power2=pow(2,num_mlfmalevels-indlev-1);
        int power8=pow(8,num_mlfmalevels-indlev);
        temp_nodes_c[0]=node_box[node_box_forsort[0].elemID*3+0]/power2;
        temp_nodes_c[1]=node_box[node_box_forsort[0].elemID*3+1]/power2;
        temp_nodes_c[2]=node_box[node_box_forsort[0].elemID*3+2]/power2;
        box_pos[(ind_filledboxes+num_filledboxes[indlev])*3+0] = boxmin[0]+temp_nodes_c[0]*size_boxes[indlev]+boxsizesvec;
        box_pos[(ind_filledboxes+num_filledboxes[indlev])*3+1] = boxmin[1]+temp_nodes_c[1]*size_boxes[indlev]+boxsizesvec;
        box_pos[(ind_filledboxes+num_filledboxes[indlev])*3+2] = boxmin[2]+temp_nodes_c[2]*size_boxes[indlev]+boxsizesvec;
        for(int indnode=1; indnode<num_nodes; indnode++){
            temp_nodes[indnode-1]=node_box_forsort[indnode-1].elemVal/power8;
            temp_nodes[indnode]=node_box_forsort[indnode].elemVal/power8;
            if (temp_nodes[1]!=temp_nodes[0]){
                if((node_box_forsort[indnode].elemVal!=node_box_forsort[indnode-1].elemVal) || (indlev==num_mlfmalevels-1))
                    num_children[ind_filledboxes+num_filledboxes[indlev]]++;
            } else{
                num_filledboxes[indlev]++;
                num_children[ind_filledboxes+num_filledboxes[indlev]]=1;
                temp_nodes_c[0]=node_box[node_box_forsort[indnode].elemID*3+0]/power2;
                temp_nodes_c[1]=node_box[node_box_forsort[indnode].elemID*3+1]/power2;
                temp_nodes_c[2]=node_box[node_box_forsort[indnode].elemID*3+2]/power2;
                box_pos[(ind_filledboxes+num_filledboxes[indlev])*3+0] = boxmin[0]+temp_nodes_c[0]*size_boxes[indlev]+boxsizesvec;
                box_pos[(ind_filledboxes+num_filledboxes[indlev])*3+1] = boxmin[1]+temp_nodes_c[1]*size_boxes[indlev]+boxsizesvec;
                box_pos[(ind_filledboxes+num_filledboxes[indlev])*3+2] = boxmin[2]+temp_nodes_c[2]*size_boxes[indlev]+boxsizesvec;
            }
        }
        ind_filledboxes+=num_filledboxes[indlev];
    }
    delete[] node_box;
    delete[] node_box_forsort;

    // Prepare the start and end index of the boxes at each level
    num_filledboxes[num_mlfmalevels*1+0]=0;
    num_filledboxes[num_mlfmalevels*2+0]=0;
    for(int indlev=1; indlev<num_mlfmalevels; indlev++){
        num_filledboxes[num_mlfmalevels*2+indlev]=num_filledboxes[num_mlfmalevels*0+indlev]+num_filledboxes[num_mlfmalevels*2+indlev-1];
        num_filledboxes[num_mlfmalevels*1+indlev]=num_filledboxes[num_mlfmalevels*2+indlev]-num_filledboxes[num_mlfmalevels*0+indlev]+1;
    }

    // Flip the arrays "num_children" and "box_pos" for MLFMA operations
    int *num_childrenTemp=new int[num_filledboxesall];
    for(int ind1=0; ind1<num_filledboxesall; ind1++){
        num_childrenTemp[ind1]=num_children[num_filledboxesall-ind1-1];
    }
    for(int indlev=0; indlev<num_mlfmalevels; indlev++){
        int ind1_start=num_filledboxes[num_mlfmalevels*1+indlev];
        int ind1_end=num_filledboxes[num_mlfmalevels*2+indlev];
        for(int ind1=0; ind1<num_filledboxes[num_mlfmalevels*0+indlev]; ind1++){
            num_children[ind1_start]=num_childrenTemp[ind1_end];
            ind1_start++;
            ind1_end--;
        }
    }
    delete[] num_childrenTemp;

    double *box_posTemp=new double[num_filledboxesall*3];
    for(int ind1=0; ind1<num_filledboxesall; ind1++){
        box_posTemp[ind1*3+0]=box_pos[(num_filledboxesall-ind1-1)*3+0];
        box_posTemp[ind1*3+1]=box_pos[(num_filledboxesall-ind1-1)*3+1];
        box_posTemp[ind1*3+2]=box_pos[(num_filledboxesall-ind1-1)*3+2];
    }
    for(int indlev=0; indlev<num_mlfmalevels; indlev++){
        int ind1_start=num_filledboxes[num_mlfmalevels*1+indlev];
        int ind1_end=num_filledboxes[num_mlfmalevels*2+indlev];
        for(int ind1=0; ind1<num_filledboxes[num_mlfmalevels*0+indlev]; ind1++){
            box_pos[ind1_start*3+0]=box_posTemp[ind1_end*3+0];
            box_pos[ind1_start*3+1]=box_posTemp[ind1_end*3+1];
            box_pos[ind1_start*3+2]=box_posTemp[ind1_end*3+2];
            ind1_start++;
            ind1_end--;
        }
    }
    delete[] box_posTemp;

    // Relate the children of the parent boxes "children" and parent of children "children2parent"
    children=new int[num_filledboxesall*2];
    children2parent=new int[num_filledboxesall];
    for(int indbox=0; indbox<num_filledboxesall; indbox++){
        children[indbox*2+0]=0;
        children[indbox*2+1]=0;
        children2parent[indbox]=-1;
    }
    for(int indlev=0; indlev<num_mlfmalevels-1; indlev++){
        int firstchild=num_filledboxes[num_mlfmalevels*1+indlev+1];
        int startbox=num_filledboxes[num_mlfmalevels*1+indlev];
        int endbox=num_filledboxes[num_mlfmalevels*2+indlev]+1;
        for(int indbox=startbox; indbox<endbox; indbox++){
            children[indbox*2+0]=firstchild;
            children[indbox*2+1]=firstchild+num_children[indbox]-1;
            int startparent=firstchild;
            int endparent=firstchild+num_children[indbox];
            for(int indparent=startparent; indparent<endparent; indparent++){
                children2parent[indparent]=indbox;
            }
            firstchild+=num_children[indbox];
        }
    }
    int firstchild=0;
    int startbox=num_filledboxes[num_mlfmalevels*2-1];
    int endbox=num_filledboxes[num_mlfmalevels*3-1]+1;
    for(int indbox=startbox; indbox<endbox; indbox++){
        children[indbox*2+0]=firstchild;
        children[indbox*2+1]=firstchild+num_children[indbox]-1;
        firstchild+=num_children[indbox];
    }


    printf("\n");
    for(int indlev=0;indlev<num_mlfmalevels;indlev++){
        int num_box=num_filledboxes[indlev];
        printf("Level: %d, Number of Boxes: %d \n", indlev, num_box);
    }
    printf("\n");


}

long int dec2bin2dec(int posx, int posy, int posz, int num_lev){
    long int power2=1;
    long int rwg_pos=0;
    for(int ind1=0;ind1<num_lev;ind1++){
        rwg_pos+=long(posz%2)*power2;
        posz/=2;
        power2*=2;
        rwg_pos+=long(posy%2)*power2;
        posy/=2;
        power2*=2;
        rwg_pos+=long(posx%2)*power2;
        posx/=2;
        power2*=2;

    }
    return rwg_pos;
}

int compare (const void *a, const void *b)
{
//    return ( (*(struct sortedElements*)a).elemVal - (*(struct sortedElements*)b).elemVal );
    if(( (*(struct sortedElements*)a).elemVal - (*(struct sortedElements*)b).elemVal )<0)
        return -1;
    if(( (*(struct sortedElements*)a).elemVal - (*(struct sortedElements*)b).elemVal )>0)
        return 1;
    return 0;

}


