/************
 *  main.c  *
 ************/


/*=================================================*
 *  Hauptprogramm zum Wavelet-Galerkin-Verfahren.  *
 *=================================================*/


#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "constants.h"
#include "intvector_pwl.h"
#include "vector2.h"
#include "vector3.h"
#include "sparse2.h"
#include "sparse.h"
#include "basis_pwl.h"
#include "dwt_pwl.h"
#include "WEM_pwl.h"
#include "WEMRHS_pwl.h"
#include "WEMPCG_pwl.h"
#include "WEMPGMRES_pwl.h"
#include "compression_pwl.h"
#include "interpolate_pwl.h"
#include "read_points.h"
#include "postproc_pwl.h"
#include "topology_pwl.h"
#include "precond_pwl.h"
#include "energy_pwl.h"
#include "volume.h"
#include "kern.h"
#include "molecule.h"


#if !defined pi
#define pi 3.1415926535897932385
#endif


int main(int argc, char* argv[])
{
    sparse2 S_i, S_e;           /* komprimierte Steifigkeitsmatrix    */
    sparse G;                   /* Massenmatrix                       */
    vector3 *P;                 /* Punkteliste der Einskalenbasis     */
    unsigned int **F;           /* Elementliste der Einskalenbasis    */
    vector3 ****T;              /* Oberflaecheninterpolation          */
    vector3 ***U;               /* Knotenpunkte                       */
    unsigned int nf;            /* Laenge von F                       */
    unsigned int np;            /* Laenge von P                       */
    unsigned int p;             /* Anzahl der Patches                 */
    unsigned int M;             /* 2^M*2^M Elemente pro Patch         */
    element_pwl *E;             /* hierarchische Elementliste         */
    wavelet_pwl *W;             /* Liste der Wavelets                 */
    double *u, *v;              /* approximierter Dichtevektor        */
    double *rhs;                /* rechte Seite des Gleichungssystems */
    double res;                 /* berechnete Austauschenergien       */
    time_t t1, t2, t3;          /* Sartzeit/Zwischenzeit/Endzeit      */
    unsigned int i, j;          /* Laufindizes                        */
    double det;                 /* Determinante im anistropen Fall    */

    char* infile = "molecule.inp";
    unsigned int CASE = 2;      /* FLAG FOR CHOICE OF BIE */

    switch (argc) {
    case 3:
        CASE = atoi(argv[2]);
    case 2:
        infile = argv[1];
    case 1:
        break;
    default :
        printf("Usage:\n"); 
        printf("pwl.x [input_file[case]]\n");
        printf("Cases\n");
        printf("1: Second kind integral equation (field)\n");
        printf("2: First kind integral equation (potential)\n");
        printf("3: Full Second kind integral equation (isotropic)\n");
        printf("4: Full Second kind integral equation (anisotropic)\n");
        exit(-1);
    }

    /* Ausgabe */
    switch (CASE) {
    case 1:
        printf("Pure Poisson Equation with 2nd kind integral equation:\n");
        printf("======================================================\n\n");
        break;
    case 2:
        printf("Pure Poisson Equation with 1st kind integral equation:\n");
        printf("======================================================\n\n");
        break;
    case 3:
        printf("Poisson-Boltzmann Equation:\n");
        printf("===========================\n\n");
        break;
    case 4:
        printf("Anisotropic Dielectric:\n");
        printf("=======================\n\n");
        break;
    default:
        printf("ERROR: non-existent scheme\n");
        return (0);
    }

    /* Initialisierung */
    read_points(&U, &p, &M);
    read_molecule("molecule.inp");
    nf = p * (1 << M) * (1 << M);       /* Anzahl der Patches */
    time(&t1);

    /* Topologie bestimmen */
    init_interpolate_pwl(&T, U, p, M);
    np = gennet_pwl(&P, &F, U, p, M);
    free_points(&U, p, M);
    volume(F, T, p, M);

    /* erstelle Element-/Waveletliste */
    printf("Number of levels:                %d \n", M);
    printf("Number of parameter domains:     %d \n", p);
    printf("Number of ansatz functions:      %d \n", np);
    printf("Computing the overhead:          ");
    generate_elementlist_pwl(&E, P, F, p, M);
    generate_waveletlist_pwl(&W, E, p, M, np);
    set_quadrature_level_pwl(W, E, p, M, np);
    simplify_waveletlist_pwl(W, E, p, M, np);
    complete_elementlist_pwl(W, E, p, M, np);
    time(&t2);
    printf("Computation time:                %g secs.\n\n", difftime(t2, t1));

    /* erstes Paar Steifigkeitsmatrizen aufstellen */
    printf("Computing the 1st pair of system matrices: \n");
    compression_pwl(&S_i, W, E, p, M, np);
    if (CASE < 3)
        WEM_pwl(&S_i, W, P, E, T, p, M, SingleLayerInt, DoubleLayerInt, 2 * pi * (1 + epsilon) / (1 - epsilon));
    else
        WEM_pwl(&S_i, W, P, E, T, p, M, SingleLayerInt, DoubleLayerInt, 2 * pi);
	fprint_sparse2(&S_i, "old.dat");
    postproc_pwl(&S_i, W, E, p, M);
    time(&t3);                  /* Zwischenzeit */
    printf("Computation time:                %g secs.\n\n", difftime(t3, t2));

    /* zweites Paar Steifigkeitsmatrizen aufstellen */
    if ((CASE == 3) || (CASE == 4)) {
        time(&t2);
        printf("Computing the 2nd pair of system matrices: \n");
        compression_pwl(&S_e, W, E, p, M, np);
        if (CASE == 3) {
            WEM_pwl(&S_e, W, P, E, T, p, M, SingleLayerExt, DoubleLayerExt, -2 * pi);
            for (i = 0; i < np; i++) {  /* correct scaling */
                for (j = 0; j < S_e.row_number[i]; j++)
                    S_e.value1[i][j] /= epsilon;
            }
        } else {
            det = sqrt(epsilon11 * epsilon22 * epsilon33 + 
                       epsilon12 * epsilon23 * epsilon31 + 
                       epsilon13 * epsilon21 * epsilon32 - 
                       epsilon11 * epsilon32 * epsilon23 - 
                       epsilon21 * epsilon12 * epsilon33 - 
                       epsilon31 * epsilon22 * epsilon13);
            WEM_pwl(&S_e, W, P, E, T, p, M, SingleLayerAni, DoubleLayerAni, -2 * pi / det);
            for (i = 0; i < np; i++) {  /* correct scaling */
                for (j = 0; j < S_e.row_number[i]; j++) {
                    S_e.value1[i][j] *= det;
                    S_e.value2[i][j] *= det;
                }
            }
        }
        postproc_pwl(&S_e, W, E, p, M);
        time(&t3);              /* Zwischenzeit */
        printf("Computation time:                %g secs.\n\n", difftime(t3, t2));
    }

    /* loese Gleichungssystem */
    u = (double *) calloc(np, sizeof(double));
    v = (double *) calloc(np, sizeof(double));
    if (CASE == 1) {
        WEMRHS1_pwl(&rhs, W, E, T, p, M, np);
        i = WEMPGMRES1_pwl(&S_i, rhs, u, eps, W, F, p, M);
        printf("Solving the linear system:       %d iterations\n", i);
    } else if (CASE == 2) {
        WEMRHS2_pwl(&rhs, W, E, T, p, M, np); //compute N_\rho  /* compute correct rhs: b-G*A2^(-1)*b */
        i = WEMPGMRES2_pwl(&S_i, rhs, v, eps, W, F, p, M); // result in v=A2^(-1)*rhs
        fprint_vec(v, np, "v_old.dat");
        printf("Solving the 1st linear system:   %d iterations\n", i);
        init_sparse(&G, np, np, 10);
        single_scale_gram_pwl(&G, F, p, M); // compute mass matrix
        tdwtLin(v, F, M, p, np); // wavelet transform
        for (i = 0; i < np; i++) {
            for (j = 0; j < G.row_number[i]; j++) {
                u[i] += G.value[i][j] * v[G.index[i][j]];
            }
        }
        dwtLin(u, F, M, p, np); //inverse transform
        for (i = 0; i < np; i++)
            rhs[i] += 4 * pi * u[i] / (epsilon - 1);  // RHS equation (2.7)
        memset(u, 0, np * sizeof(double));
        fprint_vec(rhs, np, "rhs_old2.dat");
        i = WEMPCG_pwl(&S_i, rhs, u, eps, W, F, p, M); //V \sigma = RHS
        printf("Solving the 2nd linear system:   %d iterations\n", i);
    } else if ((CASE == 3) || (CASE == 4)) {
        WEMRHS2_pwl(&rhs, W, E, T, p, M, np); // comp potential
        i = WEMPCG_pwl(&S_i, rhs, u, eps, W, F, p, M);  /* u = V_i^(-1)*N_f */
        printf("Solving the 1st linear system:   %d iterations\n", i);
        memset(rhs, 0, np * sizeof(double));
        for (i = 0; i < np; i++) {      /* rhs = V_e*u */
            for (j = 0; j < S_e.row_number[i]; j++)
                rhs[i] += S_e.value1[i][j] * u[S_e.index[i][j]];
        }
        i = WEMPGMRES3_pwl(&S_i, &S_e, rhs, v, eps, W, F, p, M);        /* solve complicated_system u = A^(-1)*rhs */
        printf("Solving the 2nd linear system:   %d iterations\n", i);
        for (i = 0; i < np; i++)
            u[i] -= 4 * pi * v[i];      /* u = u - 4*pi*v */
    }

    time(&t2);
    printf("Computation time:                %g secs.\n\n", difftime(t2, t3));

    /* Energie-Berechnung */
    tdwtLin(u, F, M, p, np);
    fprint_vec(u, np, "u_old.dat");
    res = energy_orig_pwl(u, F, T, p, M);
    time(&t2);
    printf("Over-all computation time:       %g secs.\n", difftime(t2, t1));

    /* Speicherplatz freigeben */
    printf("\n\n\n");
    if ((CASE == 3) || (CASE == 4))
        free_sparse2(&S_e);
    free_waveletlist_pwl(&W, np);
    free_elementlist_pwl(&E, p, M);
    free_interpolate_pwl(&T, p, M);
    free_patchlist_pwl(&F, nf);
    free_sparse2(&S_i);
    free(rhs);
    free(P);
    free(u);
    free(v);
    free_molecule();
    return (0);
}

