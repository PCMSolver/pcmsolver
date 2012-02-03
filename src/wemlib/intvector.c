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

signed int 
search_integral(I, i)
	intvector      *I;
	unsigned int    i;
{
	unsigned int    low, mid, high;

	if (I->integral_number == 0)
		return (-1);	/* Liste leer */

/* binary search */
	low = 0;
	high = I->integral_number - 1;
	while (low < high) {
		mid = (low + high) / 2;
		if (I->index[mid] < i)
			low = mid + 1;
		else if (I->index[mid] > i)
			high = mid;
		else
			return (mid);
	}

	if (I->index[low] == i)
		return (low);
	else
		return (-1);	/* Eintrag nicht vorhanden */
}


/*=============*
 *  I(i) := z  *
 *=============*/

void 
set_integral(I, i, z)
/* der Eintrag darf nicht vorhanden sein */
	intvector      *I;
	unsigned int    i;
	double         *z;
{
	unsigned int    rn, low, mid, high;

/* muss die Anzahl der Elemente erhoeht werden ? */
	rn = I->integral_number++;
	if (rn % delta == 0) {
		if (rn == 0) {
			I->index = (unsigned int *) malloc(delta * sizeof(unsigned int));
			I->value = (integral *) malloc(delta * sizeof(integral));
			I->index[0] = i;
			I->value[0].sub[0] = z[0];
			I->value[0].sub[1] = z[1];
			I->value[0].sub[2] = z[2];
			return;
		} else {
			I->index = (unsigned int *) realloc(I->index, (rn + delta) * sizeof(unsigned int));
			I->value = (integral *) realloc(I->value, (rn + delta) * sizeof(integral));
		}
	}
/* darf Eintrag einfach angehaengt werden ? */
	if (I->index[rn - 1] < i) {
		I->index[rn] = i;
		I->value[rn].sub[0] = z[0];
		I->value[rn].sub[1] = z[1];
		I->value[rn].sub[2] = z[2];
		return;
	}
/* suche Eintrag */
	low = 0;
	high = rn - 1;
	while (low < high) {
		mid = (low + high) / 2;
		if (I->index[mid] < i)
			low = mid + 1;
		else
			high = mid;
	}

/* schreibe Eintrag in die Liste */
	memmove(&I->index[low + 1], &I->index[low], (rn - low) * sizeof(unsigned int));
	memmove(&I->value[low + 1], &I->value[low], (rn - low) * sizeof(integral));
	I->index[low] = i;
	I->value[low].sub[0] = z[0];
	I->value[low].sub[1] = z[1];
	I->value[low].sub[2] = z[2];
	return;
}
