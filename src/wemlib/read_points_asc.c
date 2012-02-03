/***********************
 *  read_points_asc.c  *
 ***********************/


/*==============================*
 *  Liest die Punkteliste ein.  *
 *==============================*/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "vector3.h"

int 
read_points1(P, p, m, fname)
	vector3     ****P;	/* Knotenpunkte            	  	  */
	unsigned int   *p;	/* Anzahl der Patches       	  	  */
	unsigned int   *m;	/* Zahl der Level */
	char           *fname;	/* Filename of the cavity definition */
{
	unsigned int    n = 0;	/* n*n Elemente pro Patch auf dem Level m */
	unsigned int    k, l;	/* Laufindizes                            */
	unsigned int    j1, j2, j3;	/* Punktindizes                           */
	float           x, y, z;/* Punktkoordinaten 			  */
	FILE           *file;	/* Ausgabe-File                           */

/* file with geometric data */
	file = fopen(fname, "r");
	if (file == NULL) {
		return -1;
	} else {		/* Speicherplatz allokieren */
		fscanf(file, "%d\n%d\n", m, p);
		(*P) = (vector3 ***) malloc(*p * sizeof(vector3 **));
		for (k = 0; k < *p; k++) {
			n = 1 << (*m);
			(*P)[k] = (vector3 **) malloc((n + 1) * sizeof(vector3 *));
			for (l = 0; l <= n; l++)
				(*P)[k][l] = (vector3 *) malloc((n + 1) * sizeof(vector3));
		}

		/* Daten einlesen */
		for (k = 0; k < *p * (n + 1) * (n + 1); k++) {
			fscanf(file, "%d %d %d %g %g %g\n", &j1, &j3, &j2, &x, &y, &z);
			/* Punkt entspricht dem Tupel (patch,y-index,x-index) */
			(*P)[j1][j2][j3] = vector3_make(x, y, z);
		}
		fclose(file);
	}
	return 0;
}

void 
read_points(P, p, m)
	vector3     ****P;	/* Knotenpunkte            	  	  */
	unsigned int   *p;	/* Anzahl der Patches       	  	  */
	unsigned int   *m;	/* Zahl der Level */
{
	int             ret = read_points1(P, p, m, "benzene2.dat");
	return;
}

void 
free_points(P, p, m)
	vector3     ****P;	/* Knotenpunkte            	  	  */
	unsigned int    p;	/* Anzahl der Patches       	  	  */
	unsigned int    m;	/* Zahl der Level                         */
{
	unsigned int    n = 1 << m;	/* n*n Elemente pro Patch auf dem Level m */
	unsigned int    k, l;	/* Laufindizes                            */

	for (k = 0; k < p; k++) {
		for (l = 0; l <= n; l++) {
			free((*P)[k][l]);
		}
		free((*P)[k]);
	}
	free(*P);
	return;
}

void 
alloc_points(P, p, m)
	vector3     ****P;	/* Knotenpunkte            	  	  */
	unsigned int    p;	/* Anzahl der Patches       	  	  */
	unsigned int    m;	/* Zahl der Level                         */
{
	unsigned int    n = 1 << m;	/* n*n Elemente pro Patch auf dem Level m */
	unsigned int    k, l;	/* Laufindizes                            */

	(*P) = (vector3 ***) malloc(p * sizeof(vector3 **));
	for (k = 0; k < p; k++) {
		n = 1 << (m);
		(*P)[k] = (vector3 **) malloc((n + 1) * sizeof(vector3 *));
		for (l = 0; l <= n; l++) {
			(*P)[k][l] = (vector3 *) malloc((n + 1) * sizeof(vector3));
		}
	}
	return;
}
