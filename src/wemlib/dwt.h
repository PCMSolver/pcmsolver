#ifndef DWT
#define DWT
/************
 *  dwt1.h  *
 ************/
 
 
/*===========================================*
 *  Dieses Modul enthaelt alle Routinen der  *
 *  schnellen Wavelet Transformationen.      *
 *===========================================*/
 

void dwtKon(double *a, unsigned int M, unsigned int nf);
/* Diskrete Wavelet-Transformation */


void tdwtKon(double *a, unsigned int M, unsigned int nf);
/* transformierte Diskrete Wavelet-Transformation */
#endif
