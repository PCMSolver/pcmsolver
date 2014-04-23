/*==================*
 *  Mask24_mkIDK.c  *
 *==================*/


/*======================================================*
 *  Definiert die FWT-Masken fuer stueckweise lineare   *
 *  Wavelets mit 4 verschwinden Momenten und doppelten  * 
 *  Knoten an den Patchraendern.                        *
 *======================================================*/


#include <stdlib.h>
#include "sparse.h"
#include "mask_pwl.h"


const unsigned int td_pwl = 4;
const unsigned int minLevel_pwl = 2;


void mask_T22(T, m, M)
/* definiert fuer die Wavelets mit 2 verschwindenden 
   Momente die Maske T auf dem Level m */
sparse *T;
unsigned int m, M;
{
    unsigned int n = 1 << m;
    unsigned int i;

    init_sparse(T, n + 1, n + 1, 4);

    for (i = 0; i <= n; i++) {
        if (i % 2 == 0) {       /* scaling functions */
            if (i == 0) {       /* left scaling function */
                set_sparse(T, 0, 0, +1.0);
                set_sparse(T, 0, 1, +0.5);
            } else if (i < n) { /* stationary scaling function */
                set_sparse(T, i, i - 1, 0.5);
                set_sparse(T, i, i, 1.0);
                set_sparse(T, i, i + 1, 0.5);
            } else if (i == n) {        /* right scaling function */
                set_sparse(T, n, n - 1, +0.5);
                set_sparse(T, n, n, +1.0);
            }
        } else {                /* wavelet */
            if (i == 1) {       /* left boundary wavelet */
                set_sparse(T, 1, 0, -0.7500);
                set_sparse(T, 1, 1, +0.5625);
                set_sparse(T, 1, 2, -0.1250);
                set_sparse(T, 1, 3, -0.0625);
            } else if (i < n - 1) {     /* stationary wavelet */
                set_sparse(T, i, i - 1, -0.5);
                set_sparse(T, i, i, +1.0);
                set_sparse(T, i, i + 1, -0.5);
            } else if (i == n - 1) {    /* right boundary wavelet */
                set_sparse(T, n - 1, n - 3, -0.0625);
                set_sparse(T, n - 1, n - 2, -0.1250);
                set_sparse(T, n - 1, n - 1, +0.5625);
                set_sparse(T, n - 1, n, -0.7500);
            }
        }
    }
    return;
}


void mask_T24(T, m, M)
/* definiert fuer die Wavelets mit 4 verschwindenden 
   Momente die Maske T auf dem Level m */
sparse *T;
unsigned int m, M;
{
    unsigned int n = 1 << m;
    unsigned int i;

    init_sparse(T, n + 1, n + 1, 8);

    for (i = 0; i <= n; i++) {
        if (i % 2 == 0) {       /* scaling functions */
            if (i == 0) {       /* left scaling function */
                set_sparse(T, 0, 0, +1.0);
                set_sparse(T, 0, 1, +0.5);
            } else if (i < n) { /* stationary scaling function */
                set_sparse(T, i, i - 1, 0.5);
                set_sparse(T, i, i, 1.0);
                set_sparse(T, i, i + 1, 0.5);
            } else if (i == n) {        /* right scaling function */
                set_sparse(T, n, n - 1, +0.5);
                set_sparse(T, n, n, +1.0);
            }
        } else {                /* wavelet */
            if (i == 1) {       /* left boundary wavelet */
                set_sparse(T, 1, 0, -35. / 64);
                set_sparse(T, 1, 1, 875. / 1536);
                set_sparse(T, 1, 2, -241. / 768);
                set_sparse(T, 1, 3, -53. / 512);
                set_sparse(T, 1, 4, 41. / 384);
                set_sparse(T, 1, 5, 67. / 1536);
                set_sparse(T, 1, 6, -5. / 256);
                set_sparse(T, 1, 7, -5. / 512);
            } else if (i < n - 1) {     /* stationary wavelet */
                set_sparse(T, i, i - 2, +0.125);
                set_sparse(T, i, i - 1, -0.500);
                set_sparse(T, i, i, +0.750);
                set_sparse(T, i, i + 1, -0.500);
                set_sparse(T, i, i + 2, +0.125);
            } else if (i == n - 1) {    /* right boundary wavelet */
                set_sparse(T, n - 1, n - 7, -5. / 512);
                set_sparse(T, n - 1, n - 6, -5. / 256);
                set_sparse(T, n - 1, n - 5, 67. / 1536);
                set_sparse(T, n - 1, n - 4, 41. / 384);
                set_sparse(T, n - 1, n - 3, -53. / 512);
                set_sparse(T, n - 1, n - 2, -241. / 768);
                set_sparse(T, n - 1, n - 1, 875. / 1536);
                set_sparse(T, n - 1, n, -35. / 64);
            }
        }
    }
    return;
}


void dwt_mask_pwl(T, L, m, M)
/* waehlt in Abhaengigkeit von m und M die richtige Masken */
sparse *T, *L;
unsigned int m, M;
{
/* richtige Maske auswaehlen */
    if (m <= 2)
        mask_T22(T, m, M);
    else
        mask_T24(T, m, M);
    return;
}
