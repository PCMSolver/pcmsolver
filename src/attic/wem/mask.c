/*================*
 *  Mask13_mkI.c  *
 *================*/


/*=======================================================*
 *  Definiert die FWT-Masken fuer stueckweise konstante  *
 *  Tensorprodukt-Wavelets mit 3 verschwinden Momente.   *
 *=======================================================*/


#include <stdlib.h>
#include "sparse.h"
#include "mask.h"


const unsigned int td = 3;
const unsigned int minLevel = 1;


void mask_T11(T, m)
/* definiert die Maske T auf dem Level m */
sparse *T;
unsigned int m;
{
    unsigned int n = 1 << (m - 1);
    unsigned int i;

    init_sparse(T, n, 2 * n, 2);

    for (i = 0; i < n; i++) {
        set_sparse(T, i, 2 * i, +1);
        set_sparse(T, i, 2 * i + 1, -1);
    }
    return;
}


void mask_T13(T, m)
/* definiert die Maske T auf dem Level m */
sparse *T;
unsigned int m;
{
    unsigned int n = 1 << (m - 1);
    unsigned int i;

    init_sparse(T, n, 2 * n, 6);

    for (i = 0; i < n; i++) {
        if (i == 0) {           /* left boundary wavelet */
            set_sparse(T, 0, 0, +0.625);
            set_sparse(T, 0, 1, -1.375);
            set_sparse(T, 0, 2, +0.500);
            set_sparse(T, 0, 3, +0.500);
            set_sparse(T, 0, 4, -0.125);
            set_sparse(T, 0, 5, -0.125);
        } else if (i < n - 1) { /* stationary wavelet */
            set_sparse(T, i, 2 * i - 1, -0.25);
            set_sparse(T, i, 2 * i, +0.75);
            set_sparse(T, i, 2 * i + 1, -0.75);
            set_sparse(T, i, 2 * i + 2, +0.25);
        } else if (i == n - 1) {        /* right boundary wavelet */
            set_sparse(T, n - 1, 2 * n - 6, +0.125);
            set_sparse(T, n - 1, 2 * n - 5, +0.125);
            set_sparse(T, n - 1, 2 * n - 4, -0.500);
            set_sparse(T, n - 1, 2 * n - 3, -0.500);
            set_sparse(T, n - 1, 2 * n - 2, +1.375);
            set_sparse(T, n - 1, 2 * n - 1, -0.625);
        }
    }
    return;
}


void dwt_mask(T, L, m)
/* waehlt in Abhaengigkeit von m die richtige Masken */
sparse *T, *L;
unsigned int m;
{
/* richtige Maske auswaehlen */
    if (m < 3) {
        mask_T11(T, m);
    } else {
        mask_T13(T, m);
    }
    return;
}
