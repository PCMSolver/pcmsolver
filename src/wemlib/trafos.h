#ifndef TRAFOS
#define TRAFOS
/**************
 *  Trafos.h  *
 **************/


/*==================================================*
 *  Enthaelt alle Hilfsroutinen fuer die Quadratur  *
 *==================================================*/


vector2 Kappa(vector2 a, vector2 b, double h);
/* Abbildung vom Einheitsquadrat auf das Element */


vector2 Tau(double b_x, double b_y, unsigned int CASE);
/* definiert die Drehungen fuer den Duffy-Trick */


unsigned int compare(unsigned int *F1, unsigned int *F2, unsigned int *ind1, unsigned int *ind2);
/* Untersucht zwei Patches (F1[0],...,F1[3]) und (F2[0],...,F2[3])
   auf eine gemeinsame Kante (Funktionsergebnis: 3) oder einen
   gemeinsamen Punkt (Funktionsergebnis: 4), wobei in diesen Faellen
   die Drehungen ind1 und ind2 entsprechend der Situation richtig gesetzt
   werden. Im Falle keiner Gemeinsamkeiten ist das Funktionsergebnis 1. */


void quadrature_grade_smooth(signed int *g1, signed int *g2, unsigned int m1, unsigned int m2, double dist, double prec);
/* bestimmt fuer Elemente der Level m1 und m2
   mit dem Abstand dist die benoetigten Quadraturgrade
   -> beide Elemente leben auf dem gleichem Patch */


void quadrature_grade_nonsmooth(signed int *g1, signed int *g2, unsigned int m1, unsigned int m2, double dist, double prec);
/* bestimmt fuer Elemente der Level m1 und m2
   mit dem Abstand dist die benoetigten Quadraturgrade
   -> beide Elemente leben auf verschiedenen Patches */
#endif
