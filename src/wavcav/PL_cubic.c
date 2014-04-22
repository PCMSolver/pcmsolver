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

/*
 * Purpose: Patch representation of molecular cavities from
 *			atomic coordinates and radii.
 * Version: September 21, 2009.
 * Author : Maharavo Randrianarivony.
 * Affiliation: Institut für Numerische Simulation.
 *              Rheinische Friedrich-Wilhelm Universität Bonn.
 *              Wegelerstraße 6, Bonn 53115.
 *              Germany.
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cavity.h"
#include "splinemol.h"
#include "geodesic.h"
#include "pln_sph.h"
#include "sas.h"
#include "eval.h"
#include "smooth.h"


extern double prec_probe;
extern int prec_ex;


int nerv_find_wovt_nulg(double T, point * sam, int n_sam, double *len, double eps_min, double eps_max)
{
    int q = -1, nb_trials = 20, i, tr;
    double eps, diff, step, lambda;
    step = 1.0 / (double) nb_trials;
    for (tr = 0; tr <= nb_trials; tr++) {
        lambda = (double) tr *step;
        eps = lambda * eps_max + (1.0 - lambda) * eps_min;
        for (i = 0; i < n_sam - 1; i++) {
            diff = fabs(T - len[i]);
            if (diff < eps) {
                q = i;
                break;
            }
        }
        if (q != -1)
            break;
    }
    return q;
}


void vefm_eval_bilt(double t, point * sam, int n_sam, double *len, point * P)
{
    int i, q;
    double L, T, lambda, d, D;
    double eps_min = 1.0e-7, eps_max = 1.0e-3;
    if (n_sam == 1) {
        getf_find_rogc_todj(sam[0], P);
        return;
    }
    L = len[n_sam - 1];
    T = L * t;
    if (T < 0.0)
        T = 0.0;
    if (T > L)
        T = L;
    q = -1;
    for (i = 0; i < n_sam - 1; i++)
        if ((len[i] <= T) && (T <= len[i + 1])) {
            q = i;
            break;
        }
    if (q == -1)
        q = nerv_find_wovt_nulg(T, sam, n_sam, len, eps_min, eps_max);
    if (q == -1) {
        fprintf(tmpout, "Warning: unable to find parameter position\n");
        fprintf(tmpout, "T=%f\n", T);
        for (i = 0; i < n_sam; i++)
            fprintf(tmpout, "sam[%d]=[%f,%f,%f]   len[%d]=%f\n", i, sam[i].absi, sam[i].ordo, sam[i].cote, i, len[i]);
        exit(0);
    }
    d = T - len[q];
    D = len[q + 1] - len[q];
    lambda = d / D;
    P->absi = (1.0 - lambda) * sam[q].absi + lambda * sam[q + 1].absi;
    P->ordo = (1.0 - lambda) * sam[q].ordo + lambda * sam[q + 1].ordo;
    P->cote = (1.0 - lambda) * sam[q].cote + lambda * sam[q + 1].cote;
}


void saqg_cubi_wilm(point * sam, int n_sam, int n, ns_curv * B)
{
    int i;
    double *t, *len, step;
    point *X;

    len = (double *) malloc(n_sam * sizeof(double));
    len[0] = 0.0;
    for (i = 1; i < n_sam; i++)
        len[i] = len[i - 1] + wodt_dist_gilq(sam[i - 1], sam[i]);

    t = (double *) malloc((n + 1) * sizeof(double));
    X = (point *) malloc((n + 1) * sizeof(point));
    step = 1.0 / (double) n;
    for (i = 0; i <= n; i++) {
        t[i] = (double) i *step;
        vefm_eval_bilt(t[i], sam, n_sam, len, &X[i]);
    }
    free(len);

    ralc_inte_zuts(t, X, n, B);
    free(t);
    free(X);
}


void pumn_disc_hanv(ns_curv C, int N, point * p)
{
    int i;
    double a, b, step, t;
    a = C.v0;
    b = C.v1;
    step = (b - a) / ((double) N - 1.0);
    for (i = 0; i < N; i++) {
        t = a + (double) i *step;
        cuwd_eval_nivk(C, t, &p[i]);
    }
}


void fopl_simp_dufk(int N, int M, efajor_sion2D * quad)
{
    int i, j, k, **map;
    double stepx, stepy;
    map = (int **) malloc(N * sizeof(int *));
    for (i = 0; i < N; i++)
        map[i] = (int *) malloc(M * sizeof(int));

    stepx = 1.0 / ((double) N - 1.0);
    stepy = 1.0 / ((double) M - 1.0);
    k = 0;
    for (j = 0; j < M; j++)
        for (i = 0; i < N; i++) {
            quad->knot[k].u = (double) i *stepx;
            quad->knot[k].v = (double) j *stepy;
            map[i][j] = k;
            k++;
        }

    k = 0;
    for (i = 0; i < N - 1; i++)
        for (j = 0; j < M - 1; j++) {
            quad->elem[k].frvrt = map[i][j];
            quad->elem[k].scvrt = map[i + 1][j];
            quad->elem[k].thvrt = map[i + 1][j + 1];
            quad->elem[k].ftvrt = map[i][j + 1];
            k++;
        }
    quad->n_grs = N * M;
    quad->e_grs = (N - 1) * (M - 1);

    for (i = 0; i < N; i++)
        free(map[i]);
    free(map);
}


void qurf_expo_sorp(fajor_sion3D quad)
{
    int nnd, i;
    FILE *fp;
    fp = fopen("MIOTY/qud_nodes.dat", "w");
    nnd = quad.n_grs;
    for (i = 0; i < nnd; i++)
        fprintf(fp, "%f  %f   %f\n", quad.knot[i].absi, quad.knot[i].ordo, quad.knot[i].cote);
    fclose(fp);
}


void wodm_find_votw_nuwr(fajor_sion3D quad)
{
    int nel, i, n1, n2, n3, n4;
    FILE *fp;
    fp = fopen("MIOTY/qud_elements.dat", "w");
    nel = quad.e_grs;
    for (i = 0; i < nel; i++) {
        n1 = quad.elem[i].frvrt;
        n2 = quad.elem[i].scvrt;
        n3 = quad.elem[i].thvrt;
        n4 = quad.elem[i].ftvrt;
        fprintf(fp, "%d  %d  %d  %d\n", n1 + 1, n2 + 1, n3 + 1, n4 + 1);
    }
    fclose(fp);
}


void kers_find_hasn_jumv(int N, int M, int n_surf)
{
    int nel, k, i;
    double rd, gr, bl;
    FILE *fp;
    fp = fopen("MIOTY/qud_tcolor.dat", "w");
    nel = (N - 1) * (M - 1);
    for (k = 0; k < n_surf; k++) {
        tesr_colo_donr(k, &rd, &gr, &bl);
        for (i = 0; i < nel; i++)
            fprintf(fp, "%f  %f  %f\n", rd, gr, bl);
    }
    fclose(fp);
}


void talv_quad_bapt(ns_surf * surf, int n_surf, int N, int M, fajor_sion3D * QUAD)
{
    int nnd, nel, i, k, n1, n2, n3, n4, NND, NEL;
    double x, y;
    point temp;
    efajor_sion2D quad;
    quad.knot = (parm *) malloc(N * M * sizeof(parm));
    quad.elem = (efajor *) malloc((N - 1) * (M - 1) * sizeof(efajor));
    fopl_simp_dufk(N, M, &quad);
    nnd = quad.n_grs;
    nel = quad.e_grs;
    NND = 0;
    NEL = 0;
    for (k = 0; k < n_surf; k++) {
        for (i = 0; i < nel; i++) {
            n1 = quad.elem[i].frvrt;
            n2 = quad.elem[i].scvrt;
            n3 = quad.elem[i].thvrt;
            n4 = quad.elem[i].ftvrt;
            QUAD->elem[NEL].frvrt = NND + n1;
            QUAD->elem[NEL].scvrt = NND + n2;
            QUAD->elem[NEL].thvrt = NND + n3;
            QUAD->elem[NEL].ftvrt = NND + n4;
            NEL++;
        }
        for (i = 0; i < nnd; i++) {
            x = quad.knot[i].u;
            y = quad.knot[i].v;
            cilj_eval_qelf(surf[k], x, y, &temp);
            getf_find_rogc_todj(temp, &QUAD->knot[NND]);
            NND++;
        }
    }
    free(quad.elem);
    free(quad.knot);
    QUAD->n_grs = NND;
    QUAD->e_grs = NEL;
}


void mekn_expo_homd(megamanif MG)
{
    int i, j, nnd;
    FILE *fp;
    fp = fopen("MIOTY/patc_nodes.dat", "w");
    for (j = 0; j < MG.mw_grs; j++) {
        nnd = MG.msh[j].n_grs;
        for (i = 0; i < nnd; i++)
            fprintf(fp, "%f  %f  %f\n", MG.msh[j].knot[i].absi, MG.msh[j].knot[i].ordo, MG.msh[j].knot[i].cote);
    }
    fclose(fp);
}


void dokc_expo_punk(megamanif MG)
{
    int i, j, n1, n2, n3, nel, beg;
    FILE *fp;
    fp = fopen("MIOTY/patc_elements.dat", "w");
    beg = 1;
    for (j = 0; j < MG.mw_grs; j++) {
        nel = MG.msh[j].e_grs;
        for (i = 0; i < nel; i++) {
            n1 = MG.msh[j].entity[i].frvrt;
            n2 = MG.msh[j].entity[i].scvrt;
            n3 = MG.msh[j].entity[i].thvrt;
            fprintf(fp, "%d  %d  %d\n", beg + n1, beg + n2, beg + n3);
        }
        beg = beg + MG.msh[j].n_grs;
    }
    fclose(fp);
}


void cotr_expo_wuql(megamanif MG)
{
    int i, j, nel;
    double rd, gr, bl;
    FILE *fp;
    fp = fopen("MIOTY/patc_tcolor.dat", "w");
    for (j = 0; j < MG.mw_grs; j++) {
        nel = MG.msh[j].e_grs;
        tesr_colo_donr(j, &rd, &gr, &bl);
        for (i = 0; i < nel; i++)
            fprintf(fp, "%f  %f  %f\n", rd, gr, bl);
    }
    fclose(fp);
}



void sowc_vert_refd(ns_surf S, int N, int nb, PL_curve * PL)
{
    int i, j;
    double step_lrg, step_fin, u, v;
    step_lrg = 1.0 / (double) N;
    step_fin = 1.0 / ((double) nb - 1.0);
    for (i = 0; i <= N; i++) {
        u = (double) i *step_lrg;
        for (j = 0; j < nb; j++) {
            v = (double) j *step_fin;
            cilj_eval_qelf(S, u, v, &PL[i].vertex[j]);
        }
        PL[i].v_grs = nb;
    }
}


void ruhc_hori_lotq(ns_surf S, int N, int nb, PL_curve * PL)
{
    int i, j;
    double step_lrg, step_fin, u, v;
    step_lrg = 1.0 / (double) N;
    step_fin = 1.0 / ((double) nb - 1.0);
    for (i = 0; i <= N; i++) {
        v = (double) i *step_lrg;
        for (j = 0; j < nb; j++) {
            u = (double) j *step_fin;
            cilj_eval_qelf(S, u, v, &PL[i].vertex[j]);
        }
        PL[i].v_grs = nb;
    }
}


void gemq_find_pevg_cikh(PL_curve in, PL_curve * out)
{
    int i, N;
    N = in.v_grs;
    for (i = 0; i < N; i++)
        getf_find_rogc_todj(in.vertex[i], &out->vertex[i]);
    out->v_grs = N;
}


void lesm_find_vung(ns_surf * S, int n_surf, int N, int nb_fin, PL_curve * PL)
{
    int i, j, k;
    PL_curve *pl;
    pl = (PL_curve *) malloc((N + 1) * sizeof(PL_curve));
    for (j = 0; j < N + 1; j++)
        pl[j].vertex = (point *) malloc(nb_fin * sizeof(point));
    k = 0;
    for (i = 0; i < n_surf; i++) {
        fprintf(tmpout, "quilt for surf[%d / %d]\n", i, n_surf - 1);
        sowc_vert_refd(S[i], N, nb_fin, pl);
        for (j = 0; j < N + 1; j++) {
            gemq_find_pevg_cikh(pl[j], &PL[k]);
            k++;
        }
        ruhc_hori_lotq(S[i], N, nb_fin, pl);
        for (j = 0; j < N + 1; j++) {
            gemq_find_pevg_cikh(pl[j], &PL[k]);
            k++;
        }
    }
    for (j = 0; j < N + 1; j++)
        free(pl[j].vertex);
    free(pl);
}


void netv_find_pogj(ns_surf S, point * str, point * ter)
{
    double h = 1.0e-3;
    vect3D N;
    cilj_eval_qelf(S, 0.5, 0.5, str);
    gurn_norm_tegm(S, 0.5, 0.5, h, &N);
    ter->absi = str->absi + N.absi;
    ter->ordo = str->ordo + N.ordo;
    ter->cote = str->cote + N.cote;
}
/* warning-disabler-start */

#if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
#pragma GCC diagnostic pop
#elif defined(__ICC) || defined(__INTEL_COMPILER)
#pragma warning pop
#elif defined(__clang__)
#pragma clang diagnostic pop
#endif

/* warning-disabler-end */

