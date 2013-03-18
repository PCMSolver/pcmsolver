#ifndef POT
#define POT
/***********
 *  Pot.h  *
 ***********/


/*=============================================*
 * Wertet das Potential in den Punkten R aus.  *
 *=============================================*/


void pot(double **Pot, vector3 *R, unsigned int nr, double *u, vector3 ****T, unsigned int p, unsigned int m);


#endif
