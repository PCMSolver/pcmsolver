#ifndef TOPOLOGY
#define TOPOLOGY
/****************
 *  Topology.h  *
 ****************/


/*====================================*
 *  Erstellt die topologischen Daten  *
 *====================================*/


unsigned int gennet(vector3 **P, unsigned int ***F, vector3 ***T, unsigned int p, unsigned int m);
/* berechnet Punkt- und Patchliste */


void free_patchlist(unsigned int ***F, unsigned nf);
/* gibt den Speicherplatz der Patchliste frei */
#endif
