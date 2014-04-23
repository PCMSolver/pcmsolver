#ifndef POSTPROC
#define POSTPROC
/****************
 *  PostProc.h  *
 ****************/

/*===============================================*
 *  Postprocessing: setzt alle Eintraege, deren  *
 *  Betrag kleiner als ein levelabhaengiger      *
 *  Abschneideparameter ist, auf Null.           *
 *===============================================*/


double postproc(sparse2 *T, wavelet *W, element *E, unsigned int p, unsigned int M);

#endif
