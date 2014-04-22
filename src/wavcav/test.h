/* warning-disabler-start */

#if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wextra"
#elif defined(__ICC) || defined(__INTEL_COMPILER)
#pragma warning push
#pragma warning disable "-Wall"
#elif defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wall"
#pragma clang diagnostic ignored "-Weffc++"
#pragma clang diagnostic ignored "-Wextra"
#pragma clang diagnostic ignored "-Wdocumentation"
#pragma clang diagnostic ignored "-Wdeprecated-register"
#pragma clang diagnostic ignored "-Wincompatible-pointer-types"
#endif

/* warning-disabler-end */

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
/* warning-disabler-start */

#if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
#pragma GCC diagnostic pop
#elif defined(__ICC) || defined(__INTEL_COMPILER)
#pragma warning pop
#elif defined(__clang__)
#pragma clang diagnostic pop
#endif

/* warning-disabler-end */

