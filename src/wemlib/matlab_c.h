#ifndef MATLAB_C
#define MATLAB_C
/****************
 *  Matlab_C.h  *
 ****************/


/*=================================================*
 *  Schnittstelle zwischen Matlab und C:           *
 *  Konvertiert die Patch- und Punkteliste         *
 *  vom Matlab-Format ins C-Format und umgekehrt.  *
 *=================================================*/


void            c2matlab_pointlist(vector3 **P, double *PP, unsigned int np);
/* Konvertiert die vector3-Punkteliste P in eine
   (3*np,1)-REAL-Matrix PP fuer Matlab und gibt
   den Speicherplatz von P frei. */


vector3        *matlab2c_pointlist(double *PP, unsigned int np);
/* Konvertiert die (3*np,1)-REAL-Punkteliste PP
   in eine vector3-Punkteliste fuer C */


void            c2matlab_patchlist(unsigned int ***F, double *FF, unsigned int nf);
/* Konvertiert die (nf,4)-(unsigned int)-Matrix F in eine
   (4*nf,1)-REAL-Matrix FF fuer Matlab, wobei die Indizes
   um eins erhoeht werden. Der Speicherplatz von F wird
   freigegeben. */


unsigned int  **matlab2c_patchlist(double *FF, unsigned int nf);
/* Konvertiert die (4*nf,1)-REAL-Patchliste aus Matlab
   in eine (nf,4)-(unsigned int)-Matrix F fuer C,
   wobei die Indizes um 1 reduziert werden. */


void            free_patchlist(unsigned int ***F, unsigned int nf);
/* gibt den Speicherplatz der (nf,4)-(unsigned int)-Patchliste F frei */
#endif
