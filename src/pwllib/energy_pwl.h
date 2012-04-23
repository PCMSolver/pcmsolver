#ifndef ENERGY
#define ENERGY
/**************
 *  Energy.h  *
 **************/


/*=====================================*
 * Berechnet die potentielle Energie.  *
 *=====================================*/

double energy_orig_pwl(double *u, unsigned int **F, vector3 ****T, unsigned int p, unsigned int m);

double energy_pwl(double *u, double *potential, unsigned int **F, vector3 ****T, unsigned int p, unsigned int m);

double charge_pwl(double *u, double *charge, unsigned int **F, vector3 ****T, unsigned int p, unsigned int m);
#endif
