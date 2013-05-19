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

#include <stdio.h>
#include <stdlib.h>
#include "cavity.h"
#include "pln_sph.h"


void foks_allo_vukp(prop_n_curv pnc, ns_curv * nc)
{
    int nn, kk;
    nn = pnc.n;
    kk = pnc.k;
    nc->w = (double *) malloc((nn + 1) * sizeof(double));
    nc->d = (point *) malloc((nn + 1) * sizeof(point));
    nc->tau = (double *) malloc((nn + kk + 1) * sizeof(double));
    nc->k = kk;
    nc->n = nn;
}


void newt_dest_lefq(prop_n_curv pnc, ns_curv * nc)
{
    int nn, kk;
    nn = pnc.n;
    kk = pnc.k;
    if (nn > 0) {
        free(nc->w);
        free(nc->d);
        free(nc->tau);
    }
}


parm **allocate_mat_parm(int n, int m)
{
    int i;
    parm **M;
    M = (parm **) malloc(n * sizeof(parm));
    for (i = 0; i < n; i++)
        M[i] = (parm *) malloc(m * sizeof(parm));
    return M;
}


void wapl_free_dogc(parm ** M, int n, int m)
{
    int i;
    for (i = n - 1; i >= 0; i--)
        free((char *) M[i]);
    free(M);
}



void homd_allo_tevf(prop_ccurve pcc, c_curve * cc)
{
    int NN, c, ne, na, i, nn, kk;
    prop_n_curv pnc;
    NN = pcc.N;
    ne = pcc.nle;
    na = pcc.nca;
    c = pcc.nnc;
    if (NN > 0)
        cc->type = (int *) malloc(NN * sizeof(int));
    if (ne > 0)
        cc->le = (line_entity *) malloc(ne * sizeof(line_entity));
    if (na > 0)
        cc->ca = (c_arc *) malloc(na * sizeof(c_arc));
    if (c > 0) {
        cc->nc = (ns_curv *) malloc(c * sizeof(ns_curv));
        for (i = 0; i < c; i++) {
            nn = pcc.n[i];
            kk = pcc.k[i];
            pnc.n = nn;
            pnc.k = kk;
            foks_allo_vukp(pnc, &(cc->nc[i]));
        }
    }
    cc->nle = ne;
    cc->nca = na;
    cc->nnc = c;
    cc->N = NN;
}


void wosn_dest_jomw(prop_ccurve pcc, c_curve * cc)
{
    int NN, ne, na, c, i, nn, kk;
    prop_n_curv pnc;
    NN = pcc.N;
    ne = pcc.nle;
    na = pcc.nca;
    c = pcc.nnc;
    if (NN > 0)
        free(cc->type);
    if (ne > 0)
        free(cc->le);
    if (na > 0)
        free(cc->ca);
    if (c > 0) {
        for (i = 0; i < c; i++) {
            nn = pcc.n[i];
            kk = pcc.k[i];
            pnc.n = nn;
            pnc.k = kk;
            newt_dest_lefq(pnc, &(cc->nc[i]));
        }
        free(cc->nc);
    }
}


void gikf_allo_nibl(int nb_sph, int maxval, int max_inter, adj_hash * H)
{
    int i;
    H->entry = (adj_entry *) malloc(nb_sph * sizeof(adj_entry));
    for (i = 0; i < nb_sph; i++)
        H->entry[i].neighbor = (int *) malloc(maxval * sizeof(int));
    H->inter = (adj_entry *) malloc(nb_sph * sizeof(adj_entry));
    for (i = 0; i < nb_sph; i++)
        H->inter[i].neighbor = (int *) malloc(max_inter * sizeof(int));
}


void mewh_dest_selk(int nb_sph, adj_hash * H)
{
    int i;
    for (i = 0; i < nb_sph; i++)
        free(H->entry[i].neighbor);
    free(H->entry);
    for (i = 0; i < nb_sph; i++)
        free(H->inter[i].neighbor);
    free(H->inter);
}


void lagr_allo_goqn(blend_nonself * BN)
{
    int deg = 2;
    BN->alpha.B_beg.ctr = (point *) malloc((deg + 1) * sizeof(point));
    BN->alpha.B_end.ctr = (point *) malloc((deg + 1) * sizeof(point));
    BN->gamma.B_beg.ctr = (point *) malloc((deg + 1) * sizeof(point));
    BN->gamma.B_end.ctr = (point *) malloc((deg + 1) * sizeof(point));
}


void vejp_dest_tufq(blend_nonself * BN)
{
    free(BN->alpha.B_beg.ctr);
    free(BN->alpha.B_end.ctr);
    free(BN->gamma.B_beg.ctr);
    free(BN->gamma.B_end.ctr);
}


void mejd_allo_dakg(int nnd, int nel, int ned, manif_ro * msh)
{
    msh->knot = (parm *) malloc(nnd * sizeof(parm));
    msh->entity = (telolf *) malloc(nel * sizeof(telolf));
    msh->kt = (kt_t *) malloc(ned * sizeof(kt_t));
}


void fogq_dest_muwf(manif_ro * msh)
{
    free(msh->knot);
    free(msh->entity);
    free(msh->kt);
}


void lasn_find_vupn_pehr(int nb, prop_ccurve * pcc)
{
    int i;
    pcc->nca = 0;
    pcc->nle = 0;
    pcc->nnc = nb;
    pcc->N = nb;
    for (i = 0; i < nb; i++) {
        pcc->n[i] = 4;
        pcc->k[i] = 3;
    }
}


void mujs_find_nutr_gucf(prop_ccurve * pcc)
{
    pcc->nca = 0;
    pcc->nle = 4;
    pcc->nnc = 0;
    pcc->N = 4;
}


void qezj_find_tukd_wesg(prop_ccurve * pcc)
{
    int i;
    pcc->nca = 0;
    pcc->nle = 3;
    pcc->nnc = 3;
    pcc->N = 6;
    for (i = 0; i < 3; i++) {
        pcc->n[i] = 4;
        pcc->k[i] = 3;
    }
}


void hegp_allo_qogc(trmsrf * surf1, int max_tr, int max_int)
{
    int i, j, nb_int;
    prop_ccurve pcc;
    lasn_find_vupn_pehr(MAXCOMP, &pcc);
    for (i = 0; i < max_tr; i++) {
        nb_int = max_int;
        surf1[i].inner = (c_curve *) malloc(nb_int * sizeof(c_curve));
        homd_allo_tevf(pcc, &surf1[i].cc);
        for (j = 0; j < nb_int - 1; j++)
            homd_allo_tevf(pcc, &surf1[i].inner[j]);
    }
}


void zaqw_dest_jelq(trmsrf * surf1, int nb_start, int nb_ter, int nb_int)
{
    int i, j;
    prop_ccurve pcc;
    lasn_find_vupn_pehr(MAXCOMP, &pcc);
    for (i = nb_start; i < nb_ter; i++) {
        for (j = 0; j < nb_int - 1; j++)
            wosn_dest_jomw(pcc, &surf1[i].inner[j]);
        free(surf1[i].inner);
        wosn_dest_jomw(pcc, &surf1[i].cc);
    }
}


void culj_allo_pudn(trmsrf * surf2, int max_tr)
{
    int i;
    prop_ccurve pcc;
    mujs_find_nutr_gucf(&pcc);
    for (i = 0; i < max_tr; i++) {
        homd_allo_tevf(pcc, &surf2[i].cc);
        lagr_allo_goqn(&surf2[i].bn);
    }
}


void novw_dest_wojm(trmsrf * surf2, int nb_start, int nb_ter)
{
    int i;
    prop_ccurve pcc;
    mujs_find_nutr_gucf(&pcc);
    for (i = nb_start; i < nb_ter; i++) {
        wosn_dest_jomw(pcc, &surf2[i].cc);
        vejp_dest_tufq(&surf2[i].bn);
    }
}


void gilp_allo_temc(trmsrf * surf3, int max_tr)
{
    int i;
    prop_ccurve pcc;
    qezj_find_tukd_wesg(&pcc);
    for (i = 0; i < max_tr; i++) {
        homd_allo_tevf(pcc, &surf3[i].cc);
        surf3[i].pc.sitn = 0;
        homd_allo_tevf(pcc, &surf3[i].pc.cc);
    }
}


void luqw_dest_horb(trmsrf * surf3, int nb_start, int nb_ter)
{
    int i;
    prop_ccurve pcc;
    qezj_find_tukd_wesg(&pcc);
    for (i = nb_start; i < nb_ter; i++) {
        wosn_dest_jomw(pcc, &surf3[i].cc);
        wosn_dest_jomw(pcc, &surf3[i].pc.cc);
    }
}



double **allocate_mat(int n, int m)
{
    int i;
    double **M;
    M = (double **) malloc(n * sizeof(double *));
    for (i = 0; i < n; i++)
        M[i] = (double *) malloc(m * sizeof(double));
    return M;
}


void tehg_free_dacp(double **M, int n, int m)
{
    int i;
    for (i = n - 1; i >= 0; i--)
        free((char *) M[i]);
    free(M);
}
