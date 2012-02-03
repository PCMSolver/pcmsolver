/****************
 *  Topology.c  *
 ****************/


/*====================================*
 *  Erstellt die topologischen Daten  *
 *====================================*/


#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "vector2.h"
#include "vector3.h"
#include "interpolate.h"
#include "topology.h"


unsigned int    search_point(vector3 x, vector3 *P, unsigned int *pz);


unsigned int 
search_point(x, P, pz)
/* Suchfunktion zu gennet: Falls der Punkt x in
   der Punkteliste P vorhanden ist, liefert search_point den
   Index dieses Punktes, ansonsten wird x zur Liste hinzugefuegt
   und der entsprechende Index zurueckgegeben. */
	vector3         x, *P;
	unsigned int   *pz;
{
	unsigned int    k;

/* Suche in Punkteliste nach x */
	for (k = 0; (k < *pz) && (vector3_norm(vector3_sub(P[k], x)) > 1e-6); k++);

	if (k == *pz)
		P[(*pz)++] = x;	/* Punkt nicht in Punkteliste -> anhaengen */

	return (k);		/* Index des Punktes x als Funktionsergebnis */
}


unsigned int 
gennet(P, F, T, p, m)
/* Erstellt die Punkte- und Patchliste (im C-Format!)
   und liefert als Funktionsergebnis die Laenge der
   Punkteliste, die vom Geschlecht der Oberflaeche abhaengig ist. */
	vector3       **P;	/* Zeiger auf die Punkteliste 	       */
	unsigned int ***F;	/* Zeiger auf die Patchliste	       */
	vector3      ***T;	/* Knotenpunkte                        */
	unsigned int    p;	/* Anzahl der Patches                  */
	unsigned int    m;	/* 2^m*2^m Patches pro Parametergebiet */
{
	unsigned int    n = 1 << m;	/* n = 2^m			       */
	unsigned int    i1, i2, i3;	/* Laufindizes 			       */
	unsigned int    pz, fz, kz;	/* Punkte-/Patch-/Kantenzaehler        */
	unsigned int  **K;	/* Kantenliste                         */

/* Speicherplatz allokieren: worst case */
	K = (unsigned int **) malloc(4 * p * sizeof(unsigned int *));
	(*P) = (vector3 *) malloc(p * (n + 1) * (n + 1) * sizeof(vector3));
	(*F) = (unsigned int **) malloc(p * n * n * sizeof(unsigned int *));
	for (i1 = 0; i1 < p * n * n; i1++) {
		(*F)[i1] = (unsigned int *) calloc(4, sizeof(unsigned int));
	}
	/* Eckpunkte und Kanten der Parametergebiete bestimmen */
	pz = kz = 0;
	for (i1 = 0; i1 < p; i1++) {
		fz = i1 * n * n;

		/* Bestimme die 4 Eckpunkte des Patches */
		(*F)[fz][0] = search_point(T[i1][0][0], *P, &pz);
		(*F)[fz + n - 1][1] = search_point(T[i1][0][n], *P, &pz);
		(*F)[fz + n * n - 1][2] = search_point(T[i1][n][n], *P, &pz);
		(*F)[fz + n * (n - 1)][3] = search_point(T[i1][n][0], *P, &pz);

		/* Bestimme 1. Kante */
		for (i2 = 0; (i2 < kz) && (((*F)[fz][0] != K[i2][n]) || ((*F)[fz + n - 1][1] != K[i2][0])); i2++);
		if (i2 == kz) {	/* Kante noch nicht vorhanden -> unterteilen */
			K[kz] = (unsigned int *) malloc((n + 1) * sizeof(unsigned int));
			K[kz][0] = (*F)[fz][0];
			K[kz][n] = (*F)[fz + n - 1][1];
			for (i3 = 1; i3 < n; i3++) {
				(*P)[pz] = T[i1][0][i3];
				(*F)[fz + i3 - 1][1] = (*F)[fz + i3][0] = K[kz][i3] = pz++;
			}
			kz++;
		} else {	/* Kante schon vorhanden */
			for (i3 = 1; i3 < n; i3++)
				(*F)[fz - 1 + i3][1] = (*F)[fz + i3][0] = K[i2][n - i3];
		}

		/* Bestimme 2. Kante */
		for (i2 = 0; (i2 < kz) && (((*F)[fz + n - 1][1] != K[i2][n]) || ((*F)[fz + n * n - 1][2] != K[i2][0])); i2++);
		if (i2 == kz) {	/* Kante noch nicht vorhanden -> unterteilen */
			K[kz] = (unsigned int *) malloc((n + 1) * sizeof(unsigned int));
			K[kz][0] = (*F)[fz + n - 1][1];
			K[kz][n] = (*F)[fz + n * n - 1][2];
			for (i3 = 1; i3 < n; i3++) {
				(*P)[pz] = T[i1][i3][n];
				(*F)[fz + i3 * n - 1][2] = (*F)[fz + (i3 + 1) * n - 1][1] = K[kz][i3] = pz++;
			}
			kz++;
		} else {	/* Kante schon vorhanden */
			for (i3 = 1; i3 < n; i3++)
				(*F)[fz + i3 * n - 1][2] = (*F)[fz + (i3 + 1) * n - 1][1] = K[i2][n - i3];
		}

		/* Bestimme 3. Kante */
		for (i2 = 0; (i2 < kz) && (((*F)[fz + n * n - 1][2] != K[i2][n]) || ((*F)[fz + n * (n - 1)][3] != K[i2][0])); i2++);
		if (i2 == kz) {	/* Kante noch nicht vorhanden -> unterteilen */
			K[kz] = (unsigned int *) malloc((n + 1) * sizeof(unsigned int));
			K[kz][0] = (*F)[fz + n * n - 1][2];
			K[kz][n] = (*F)[fz + n * (n - 1)][3];
			for (i3 = 1; i3 < n; i3++) {
				(*P)[pz] = T[i1][n][n - i3];
				(*F)[fz + n * n - i3][3] = (*F)[fz + n * n - i3 - 1][2] = K[kz][i3] = pz++;
			}
			kz++;
		} else {	/* Kante schon vorhanden */
			for (i3 = 1; i3 < n; i3++)
				(*F)[fz + n * n - i3][3] = (*F)[fz + n * n - i3 - 1][2] = K[i2][n - i3];
		}

		/* Bestimme 4. Kante */
		for (i2 = 0; (i2 < kz) && (((*F)[fz + n * (n - 1)][3] != K[i2][n]) || ((*F)[fz][0] != K[i2][0])); i2++);
		if (i2 == kz) {	/* Kante noch nicht vorhanden -> unterteilen */
			K[kz] = (unsigned int *) malloc((n + 1) * sizeof(unsigned int));
			K[kz][0] = (*F)[fz + n * (n - 1)][3];
			K[kz][n] = (*F)[fz][0];
			for (i3 = 1; i3 < n; i3++) {
				(*P)[pz] = T[i1][n - i3][0];
				(*F)[fz + n * (n - i3 - 1)][3] = (*F)[fz + n * (n - i3)][0] = K[kz][i3] = pz++;
			}
			kz++;
		} else {	/* Kante schon vorhanden */
			for (i3 = 1; i3 < n; i3++)
				(*F)[fz + n * (n - i3 - 1)][3] = (*F)[fz + n * (n - i3)][0] = K[i2][n - i3];
		}

		/* Verfeinere das Patch */
		for (i2 = 1; i2 < n; i2++) {
			for (i3 = 1; i3 < n; i3++) {
				(*P)[pz] = T[i1][i2][i3];
				(*F)[fz + n * (i2 - 1) + i3 - 1][2] = (*F)[fz + n * (i2 - 1) + i3][3] = (*F)[fz + n * i2 + i3 - 1][1] = (*F)[fz + n * i2 + i3][0] = pz++;
			}
		}
	}

/* Speicherplatz wieder freigeben */
	//for (i1 = 0; i1 < kz; i1++)
		free(K[i1]);
	//free(K);

	return (pz);
}


void 
free_patchlist(F, nf)
/* Gibt den Speicherplatz der (nf,4)-(unsigned int)-Patchliste F frei */
	unsigned int ***F;
	unsigned int    nf;
{
	unsigned int    k;
	for (k = 0; k < nf; k++)
		free((*F)[k]);
	free(*F);
	return;
}
