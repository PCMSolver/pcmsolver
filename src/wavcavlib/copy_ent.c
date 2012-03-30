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

#include "cavity.h"
#include "pln_sph.h"


void cunl_find_qedf_rewn(parm p, parm * q)
{
    q->u = p.u;
    q->v = p.v;
}


void getf_find_rogc_todj(point p, point * q)
{
    q->absi = p.absi;
    q->ordo = p.ordo;
    q->cote = p.cote;
}


void sihr_find_kocf_gonh(line_entity L, line_entity * M)
{
    getf_find_rogc_todj(L.p1, &(M->p1));
    getf_find_rogc_todj(L.p2, &(M->p2));
}


void jufw_find_lukj_rapg(xform x_in, xform * x_out)
{
    int i, j;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            x_out->R[i][j] = x_in.R[i][j];
    for (i = 0; i < 3; i++)
        x_out->T[i] = x_in.T[i];
}


void zobm_find_wumq_kihf(ns_curv nci, ns_curv * nco)
{
    int nn, kk, i;
    nn = nci.n;
    kk = nci.k;
    nco->n = nn;
    nco->k = kk;

    for (i = 0; i <= nn; i++) {
        nco->w[i] = nci.w[i];
        getf_find_rogc_todj(nci.d[i], &(nco->d[i]));
    }

    for (i = 0; i <= nn + kk; i++)
        nco->tau[i] = nci.tau[i];

    nco->v0 = nci.v0;
    nco->v1 = nci.v1;

    nco->prop1 = nci.prop1;
    nco->prop2 = nci.prop2;
    nco->prop3 = nci.prop3;
    nco->prop4 = nci.prop4;

    getf_find_rogc_todj(nci.nrml, &(nco->nrml));
}


void zikt_find_jotz_jewb(trm_sph ts_in, trm_sph * ts_out)
{
    getf_find_rogc_todj(ts_in.zent, &ts_out->zent);
    ts_out->beta = ts_in.beta;
    ts_out->rad = ts_in.rad;
    ts_out->nrml.absi = ts_in.nrml.absi;
    ts_out->nrml.ordo = ts_in.nrml.ordo;
    ts_out->nrml.cote = ts_in.nrml.cote;
}


void kawl_find_wuqg_kuwg(c_arc cai, c_arc * cao)
{
    cao->ZT = cai.ZT;
    getf_find_rogc_todj(cai.p1, &(cao->p1));
    getf_find_rogc_todj(cai.p2, &(cao->p2));
    getf_find_rogc_todj(cai.p3, &(cao->p3));
    cao->A = cai.A;
    cao->B = cai.B;
    jufw_find_lukj_rapg(cai.TR, &cao->TR);
}


void kotg_find_wuhk_kemt(c_curve cci, c_curve * cco)
{
    int i, NN, ne, na, nc;
    NN = cci.N;
    ne = cci.nle;
    na = cci.nca;
    nc = cci.nnc;
    for (i = 0; i < NN; i++)
        cco->type[i] = cci.type[i];
    for (i = 0; i < ne; i++)
        sihr_find_kocf_gonh(cci.le[i], &(cco->le[i]));
    for (i = 0; i < na; i++)
        kawl_find_wuqg_kuwg(cci.ca[i], &(cco->ca[i]));
    for (i = 0; i < nc; i++)
        zobm_find_wumq_kihf(cci.nc[i], &(cco->nc[i]));
    cco->N = NN;
    cco->nle = ne;
    cco->nca = na;
    cco->nnc = nc;
}


void romh_find_cont_qucr(pt_tor P, pt_tor * Q)
{
    poms_find_resk_lonb(P.alpha, &Q->alpha);
    poms_find_resk_lonb(P.beta, &Q->beta);
    poms_find_resk_lonb(P.gamma, &Q->gamma);
    poms_find_resk_lonb(P.delta, &Q->delta);
}


void neqg_find_lodr_bogm(sphere S_in, sphere * S_out)
{
    S_out->zent.absi = S_in.zent.absi;
    S_out->zent.ordo = S_in.zent.ordo;
    S_out->zent.cote = S_in.zent.cote;
    S_out->rad = S_in.rad;
}
