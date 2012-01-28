/*****************
 *  IntVector.h  *
 *****************/
 

/*================================*
 *  Definiert den Integralvektor  *  
 *================================*/


#include <string.h>
#include <stdlib.h>
#include "intvector.h"
#include "constants.h"


/*===============================================*
 *  Suchalgorithmus gemaess binary-search:       *
 *  Liefert den Index des gewuenschten Elements  *
 *  bzw. -1 falls es nicht vorhanden.            *
 *===============================================*/

signed int search_integral(I,i)
intvector	*I;
unsigned int	i;
{
unsigned int	low, mid, high;

if (I->integral_number == 0) return(-1);	/* Liste leer */

/* binary search */
low = 0;
high = I->integral_number-1;
while (low < high)
{  mid = (low + high)/2;
   if      (I->index[mid] < i)  low = mid+1;
   else if (I->index[mid] > i)  high = mid;
   else return(mid);
   }

if (I->index[low] == i) return(low);
else                    return(-1); 	/* Eintrag nicht vorhanden */
}


/*=============*
 *  I(i) := z  *
 *=============*/

void set_integral(I,i,z)
/* der Eintrag darf nicht vorhanden sein */
intvector       *I;
unsigned int    i;
double		*z;
{
unsigned int	rn, low, mid, high;

/* muss die Anzahl der Elemente erhoeht werden ? */
rn = I->integral_number++;
if (rn%delta == 0)
{  if (rn == 0)
   {  I->index = (unsigned int*) malloc(delta*sizeof(unsigned int));
      I->value = (integral*) malloc(delta*sizeof(integral));
      I->index[0] = i;
      memcpy(I->value[0].sub,z,48*sizeof(double));
      return;
      }
   else
   {  I->index = (unsigned int*) realloc(I->index,(rn+delta)*sizeof(unsigned int));
      I->value = (integral*) realloc(I->value,(rn+delta)*sizeof(integral));
      }
   }

/* darf Eintrag einfach angehaengt werden ? */
if (I->index[rn-1] < i)
{  I->index[rn] = i;
   memcpy(I->value[rn].sub,z,48*sizeof(double));
   return;
   }

/* suche Eintrag */
low = 0;
high = rn-1;
while (low < high)
{  mid = (low+high)/2;
   if (I->index[mid] < i) low = mid+1;
   else                   high = mid;
   }
   
/* schreibe Eintrag in die Liste */
memmove(&I->index[low+1],&I->index[low],(rn-low)*sizeof(unsigned int));
memmove(&I->value[low+1],&I->value[low],(rn-low)*sizeof(integral));
I->index[low] = i;
memcpy(I->value[low].sub,z,48*sizeof(double));
return;
}


/*=======================================*
 *  bestimmt den transponierten Eintrag  *
 *=======================================*/

void permutate(a,b)
double		*a, *b;
{
a[ 0] = b[ 0];
a[ 1] = b[ 4];
a[ 2] = b[ 8];
a[ 3] = b[12];
a[ 4] = b[ 1];
a[ 5] = b[ 5];
a[ 6] = b[ 9];
a[ 7] = b[13];
a[ 8] = b[ 2];
a[ 9] = b[ 6];
a[10] = b[10];
a[11] = b[14];
a[12] = b[ 3];
a[13] = b[ 7];
a[14] = b[11];
a[15] = b[15];

a[16] = b[32];
a[17] = b[36];
a[18] = b[40];
a[19] = b[44];
a[20] = b[33];
a[21] = b[37];
a[22] = b[41];
a[23] = b[45];
a[24] = b[34];
a[25] = b[38];
a[26] = b[42];
a[27] = b[46];
a[28] = b[35];
a[29] = b[39];
a[30] = b[43];
a[31] = b[47];

a[32] = b[16];
a[33] = b[20];
a[34] = b[24];
a[35] = b[28];
a[36] = b[17];
a[37] = b[21];
a[38] = b[25];
a[39] = b[29];
a[40] = b[18];
a[41] = b[22];
a[42] = b[26];
a[43] = b[30];
a[44] = b[19];
a[45] = b[23];
a[46] = b[27];
a[47] = b[31];
}
