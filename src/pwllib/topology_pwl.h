#ifndef TOPOLOGY_PWL
#define TOPOLOGY_PWL
/****************
 *  Topology.h  *
 ****************/
#ifdef __cplusplus
extern "C" {
#endif


/*====================================*
 *  Erstellt die topologischen Daten  *
 *====================================*/


    unsigned int gennet_pwl(vector3 **P, unsigned int ***F, vector3 ***T, unsigned int p, unsigned int m);
/* berechnet Punkt- und Patchliste */


    void free_patchlist_pwl(unsigned int ***F, unsigned nf);
/* gibt den Speicherplatz der Patchliste frei */
#ifdef __cplusplus
}
#endif
#endif
