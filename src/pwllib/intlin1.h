#ifndef INTLIN1
#define INTLIN1
/***************
 *  IntLin1.h  *
 ***************/


/*===========================================================*
 *  Fernfeld-Quadratur-Routine:			             *
 *  Die in der Funktion init_randwerte vorab berechneten     *
 *  Auswertepunkte und Gewichte der Gauss-Quadratur werden   *	
 *  in IntLin1 zum entsprechenden Integral zusammengefuegt.  *
 *===========================================================*/


typedef struct {
    vector3 *Chi, *n_Chi_pwl;
    double *det_dChi;
    unsigned int nop;
} randwerte;


void init_randwerte(randwerte **RW, unsigned int g_max);
/* Initialisiert die Randwerte */


void reset_randwerte(randwerte *RW, unsigned int g_max);
/* Resetted die Randwerte */


void free_randwerte(randwerte **RW, unsigned int g_max);
/* Gibt den Speicherplatz fuer die Randwerte frei */


void IntLin1(double *c, element_pwl * element1, element_pwl * element2, randwerte *RW, cubature *Q1, cubature *Q2, vector3 ****P, unsigned int M, double SL(), double DL());
/* No-Problem-Quadrature-Routine */
#endif
