#ifndef CUBATURE
#define CUBATURE
/****************
 *  Cubature.h  *
 ****************/


/*==============================*
 *  Deklariert Kubatur-Formeln  *
 *  auf dem Referenzgebiet.     *
 *==============================*/


typedef struct {
    unsigned int nop;           /* Zahl der Stuetzstellen und Gewichte */
    vector2 *xi;                /* Stuetzstellen der Kubatur-Formel    */
    double *w;                  /* Gewichte der Kubatur-Formel         */
} cubature;

#endif
