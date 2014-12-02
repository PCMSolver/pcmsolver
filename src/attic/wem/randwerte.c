/***************
 *  IntKon1.c  *
 ***************/


/*===========================================================*
 *  Fernfeld-Quadratur-Routine:			             *
 *  Die in der Funktion init_randwerte vorab berechneten     *
 *  Auswertepunkte und Gewichte der Gauss-Quadratur werden   *
 *  in IntKon1 zum entsprechenden Integral zusammengefuegt.  *
 *===========================================================*/


#include <stdlib.h>
#include "vector3.h"
#include "randwerte.h"


void init_randwerte(RW, g_max)
/* Initialisiert die Randwerte */
randwerte **RW;
unsigned int g_max;
{
    unsigned int i;

    (*RW) = (randwerte *) malloc(g_max * sizeof(randwerte));
    for (i = 0; i < g_max; i++) {
        (*RW)[i].Chi = (vector3 *) malloc((i + 1) * (i + 1) * sizeof(vector3));
        (*RW)[i].n_Chi = (vector3 *) malloc((i + 1) * (i + 1) * sizeof(vector3));
        (*RW)[i].det_dChi = (double *) malloc((i + 1) * (i + 1) * sizeof(double));
    }
    return;
}


void reset_randwerte(RW, g_max)
/* Resetted die Randwerte */
randwerte *RW;
unsigned int g_max;
{
    unsigned int i;
    for (i = 0; i < g_max; i++)
        RW[i].nop = 0;
    return;
}


void free_randwerte(RW, g_max)
/* Gibt den Speicherplatz der Randwerte frei */
randwerte **RW;
unsigned int g_max;
{
    unsigned int i;

    for (i = 0; i < g_max; i++) {
        free((*RW)[i].Chi);
        free((*RW)[i].n_Chi);
        free((*RW)[i].det_dChi);
    }
    free(*RW);
    return;
}


