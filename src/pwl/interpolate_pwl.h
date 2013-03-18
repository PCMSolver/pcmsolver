#ifndef INTERPOLATE_PWL
#define INTERPOLATE_PWL
/*******************
 *  Interpolate.h  *
 *******************/


/*=============================================================*
 *  Enthaelt alle Routinen zur Interpolation der Oberflaeche.  *
 *=============================================================*/


void init_interpolate_pwl(vector3 *****P, vector3 ***Q, unsigned int p, unsigned int m);
/* berechnet die Koeffizienten des Interpolationspolynoms */


void free_interpolate_pwl(vector3 *****P, unsigned int p, unsigned int m);
/* gibt die Koeffizienten des Interpolationspolynoms frei */


vector3 Chi_pwl(vector2 a, vector3 ***p, unsigned int m);
/* wertet das durch p gegebene Interpolationspolynom in a aus */


vector3 dChi_dx_pwl(vector2 a, vector3 ***p, unsigned int m);
/* definiert die Ableitung nach x des Interpolationspolynoms p */


vector3 dChi_dy_pwl(vector2 a, vector3 ***p, unsigned int m);
/* definiert die Ableitung nach y des Interpolationspolynoms p */


vector3 n_Chi_pwl(vector2 a, vector3 ***p, unsigned int m);
/* definiert die Normalableitung des Interpolationspolynoms p */
#endif
