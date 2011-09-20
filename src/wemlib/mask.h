#ifndef MASK
#define MASK
/*==========*
 *  Mask.h  *
 *==========*/


/*================================================================*
 *  Definiert die FWT-Masken fuer stueckweise konstante Wavelets  *
 *================================================================*/
 

void dwt_mask(sparse *T, sparse *L, unsigned int m);
/* waehlt in Abhaengigkeit von m die richtige Maske */
#endif
