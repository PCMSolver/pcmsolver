#ifndef RANDWERTE_H_
#define RANDWERTE_H_
/***************
 *  randwerte.h  *
 ***************/

typedef struct {
    vector3 *Chi, *n_Chi;
    double *det_dChi;
    unsigned int nop;
} randwerte;


void init_randwerte(randwerte **RW, unsigned int g_max);
/* Initialisiert die Randwerte */


void reset_randwerte(randwerte *RW, unsigned int g_max);
/* Resetted die Randwerte */


void free_randwerte(randwerte **RW, unsigned int g_max);
/* Gibt den Speicherplatz fuer die Randwerte frei */
#endif
