#ifndef DWT_PWL
#define DWT_PWL
/************
 *  dwt1.h  *
 ************/


/*===========================================*
 *  Dieses Modul enthaelt alle Routinen der  *
 *  schnellen Wavelet Transformationen.      *
 *===========================================*/


void multiple(unsigned int ****C, unsigned int **Z, unsigned int **F, unsigned int M, unsigned int p, unsigned int np);
/* Hilfsfunktion zur Wavelet-Transformation: Berechnet aus der 
   Elementliste F eine Liste der lokalen Gitterpunkte (Basisliste) 
   und eine Vielfachheitenliste der Punkte */


void dwtLin(double *a, unsigned int **F, unsigned int M, unsigned int p, unsigned int np);
/* Diskrete Wavelet-Transformation */


void tdwtLin(double *a, unsigned int **F, unsigned int M, unsigned int p, unsigned int np);
/* transformierte Diskrete Wavelet-Transformation */
#endif
