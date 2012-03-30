/************
 *  molecule.c  *
 ************/


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "vector3.h"
#include "molecule.h"


/* number of charges */
static unsigned int nAtoms;

/* position of charges: vectors (x,y,z) */
static vector3 * center;

/* charges: double values */
static double * charge;

int read_molecule(char * filename)
{
    FILE* fp = fopen(filename, "r");
    if (fp == NULL) {
        return -1;
    } else {                    /* Speicherplatz allokieren */
        fscanf(fp, "%d\n", &nAtoms);
        center = (vector3 *) malloc(nAtoms * sizeof(vector3));
        charge = (double *) malloc(nAtoms * sizeof(double));
        for (int i = 0; i < nAtoms; i++) {
            fscanf(fp, "%lg %lg %lg %lg\n", &(charge[i]), &(center[i].x),  &(center[i].y),  &(center[i].z));
        }
    }
    fclose(fp);
    return 0;    
}

double potmol(vector3 a)
{
    double pot = 0;
    for (int i = 0; i < nAtoms; i++)
        pot += charge[i] / vector3_norm(vector3_make(a.x - center[i].x, a.y - center[i].y, a.z - center[i].z));
    return (pot);
}


vector3 field(vector3 a)
{
    vector3 c, r, v;
    c.x = c.y = c.z = 0;
    for (int i = 0; i < nAtoms; i++) {
        r = vector3_make(a.x - center[i].x, a.y - center[i].y, a.z - center[i].z);
        v = vector3_Smul(charge[i] / pow(vector3_norm(r), 3), r);
        c.x -= v.x;
        c.y -= v.y;
        c.z -= v.z;
    }
    return (c);
}

void free_molecule() {
    free(center);
    free(charge);
}
