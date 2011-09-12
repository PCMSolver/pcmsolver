/*******************
 *  Interpolate.h  *
 *******************/
#ifdef __cplusplus
extern "C" {
#endif


/*=============================================================*
 *  Enthaelt alle Routinen zur Interpolation der Oberflaeche.  *
 *=============================================================*/
 

void init_interpolate(vector3 *****P, vector3 ***Q, unsigned int p, unsigned int m);
/* berechnet die Koeffizienten des Interpolationspolynoms */


void free_interpolate(vector3 *****P, unsigned int p, unsigned int m);
/* gibt die Koeffizienten des Interpolationspolynoms frei */


vector3 Chi(vector2 a, vector3 ***p, unsigned int m);
/* wertet das durch p gegebene Interpolationspolynom in a aus */


vector3 dChi_dx(vector2 a, vector3 ***p, unsigned int m);
/* definiert die Ableitung nach x des Interpolationspolynoms p */


vector3 dChi_dy(vector2 a, vector3 ***p, unsigned int m);
/* definiert die Ableitung nach y des Interpolationspolynoms p */


vector3 n_Chi(vector2 a, vector3 ***p, unsigned int m);
/* definiert die Normalableitung des Interpolationspolynoms p */

#ifdef __cplusplus
}
#endif
