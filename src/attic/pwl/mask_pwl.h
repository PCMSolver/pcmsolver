#ifndef MASK_PWL
#define MASK_PWL
/*==========*
 *  Mask.h  *
 *==========*/


/*================================================================*
 *  Definiert die FWT-Masken fuer stueckweise konstante Wavelets  *
 *================================================================*/


void dwt_mask_pwl(sparse *T, sparse *L, unsigned int m, unsigned int M);
/* waehlt in Abhaengigkeit von m die richtige Maske */
#endif
