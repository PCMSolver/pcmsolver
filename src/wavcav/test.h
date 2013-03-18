#ifndef TEST_H_
#define TEST_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>
#include "cavity.h"
#include "pln_sph.h"
#include "sas.h"
#include "splinemol.h"
#include "meshsas.h"
#include "partcas.h"
#include "geodesic.h"
#include "smooth.h"


void bugt_clea_pujv(supp_sph_tri * treb, int nb_surf3);

void tofr_sort_mebd(atom * A, adj_hash H, int p);

void hetc_glob_kaqc(double probe, double c_param, sphere * S, int nb_sph, adj_hash H, int *forc_term);

int bihp_test_romg(char *filename, double probe, double c_param, int *n_size);

void cavity_create_(double *probe, double *coarsity, int *pl, int *info);

int waveletCavityDrv_(double probeRadius, double coarsity, int patchLevel, 
                      const char* infile);

#endif
