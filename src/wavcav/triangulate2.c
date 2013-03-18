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
#include "pln_sph.h"
#include "triang.h"


void wosr_find_nozs(parm org, parm A, parm * res)
{
    vect U;
    cuwl_unit_pist(org, A, &U);
    res->u = org.u - U.u;
    res->v = org.v - U.v;
}



int nojp_test_keqb(mult_conn P, int *S, int lS, int i, int pos, int k)
{
    int en, res, ni, pi, p, ts, j, b, e, nx, ind;
    int nx_k, pr_k, ts_enc;
    parm ni_t, pi_t;
    ni = rupj_next_qejp(P, i, pos);
    pi = murc_prec_kotq(P, i, pos);
    wosr_find_nozs(P.vertex[i], P.vertex[ni], &ni_t);
    wosr_find_nozs(P.vertex[i], P.vertex[pi], &pi_t);

    pr_k = k - 1;
    if (pr_k == -1)
        pr_k = lS - 1;
    nx_k = k + 1;
    if (nx_k == lS)
        nx_k = 0;
    ts_enc = kozf_test_luzk(P.vertex[S[k]], P.vertex[S[nx_k]], P.vertex[S[pr_k]], P.vertex[i]);
    if (ts_enc == 0)
        return 0;

    p = S[k];
    en = kozf_test_luzk(P.vertex[i], pi_t, ni_t, P.vertex[p]);
    if (en == 0)
        res = 0;
    else {
        ind = 1;
        for (j = 0; j < lS; j++) {
            if (j == lS - 1)
                nx = 0;
            else
                nx = j + 1;
            b = S[j];
            e = S[nx];
            if ((p != b) && (p != e) && (i != b) && (i != e)) {
                ts = pumj_test_rewf(P.vertex[p], P.vertex[i], P.vertex[b], P.vertex[e]);
                if (ts == 1) {
                    ind = 2;
                    break;
                }
            }
        }
        if (ind == 1)
            res = 1;
        else
            res = 0;
    }
    return res;
}



int qemv_test_nekr(mult_conn P, int *S, int lS, int i, int pos, int k)
{
    int en, res, ni, pi, p, ts, j, b, e, nx, ind;
    ni = rupj_next_qejp(P, i, pos);
    pi = murc_prec_kotq(P, i, pos);
    p = S[k];
    en = kozf_test_luzk(P.vertex[i], P.vertex[ni], P.vertex[pi], P.vertex[p]);
    if (en == 0)
        res = 0;
    else {
        ind = 1;
        for (j = 0; j < lS; j++) {
            if (j == lS - 1)
                nx = 0;
            else
                nx = j + 1;
            if ((j != k) && (nx != k)) {
                b = S[j];
                e = S[nx];
                ts = pumj_test_rewf(P.vertex[p], P.vertex[i], P.vertex[b], P.vertex[e]);
                if (ts == 1) {
                    ind = 2;
                    break;
                }
            }
        }
        if (ind == 1)
            res = 1;
        else
            res = 0;
    }
    return res;
}



double jegd_find_kuhw(mult_conn P, int *S, int lengthS, int j, int k, int pos)
{
    int nx, pr, pj, nj, i;
    double theta, *alpha;
    alpha = (double *) malloc(8 * sizeof(double));
    if (k == lengthS - 1)
        nx = S[0];
    else
        nx = S[k + 1];
    if (k == 0)
        pr = S[lengthS - 1];
    else
        pr = S[k - 1];
    nj = rupj_next_qejp(P, j, pos);
    pj = murc_prec_kotq(P, j, pos);
    alpha[0] = lomn_inte_cubq(P.vertex[pj], P.vertex[j], P.vertex[S[k]]);
    alpha[1] = lomn_inte_cubq(P.vertex[S[k]], P.vertex[j], P.vertex[nj]);
    alpha[2] = lomn_inte_cubq(P.vertex[pr], P.vertex[S[k]], P.vertex[j]);
    alpha[3] = lomn_inte_cubq(P.vertex[j], P.vertex[S[k]], P.vertex[nx]);
    for (i = 4; i < 8; i++)
        alpha[i] = fabs(MY_PI - alpha[i - 4]);
    theta = 10.0;
    for (i = 0; i < 8; i++)
        if (theta > alpha[i])
            theta = alpha[i];
    free(alpha);
    return theta;
}



int zotj_find_vemk(mult_conn P, int *S, int lS, int f, int j, int *l, double *angle, double beta, int search_type)
{
    int res, ind, ts, i, d, k;
    double lg, theta;
    lg = -1.0e+4;
    ind = 1;
    for (i = 0; i < lS; i++) {
        ts = 0;
        if (search_type == RESTRICTED)
            d = nojp_test_keqb(P, S, lS, j, f, i);
        if (search_type == LIBERAL)
            d = qemv_test_nekr(P, S, lS, j, f, i);
        if (d == 1) {
            theta = jegd_find_kuhw(P, S, lS, j, i, f);
            if (theta > lg) {
                lg = theta;
                k = S[i];
                *l = i;
                ind = 2;
                *angle = lg;
            }
        }
    }

    if (ind == 1)
        res = FAILURE;
    else {
        if ((lg < beta) && (search_type == RESTRICTED))
            res = FAILURE;
        else
            res = SUCCESS;
    }
    return res;
}



void kupq_find_vifr(mult_conn P, int *S, int lS, int f, int j, int *l, double *angle, double beta)
{
    int suc;
    suc = zotj_find_vemk(P, S, lS, f, j, l, angle, beta, RESTRICTED);
    if (suc == FAILURE)
        zotj_find_vemk(P, S, lS, f, j, l, angle, beta, LIBERAL);
}



int legh_find_qijk(mult_conn P, int *remain, int lg_rem, int type, int *S, int lS, double *ratio, double *angle, int *L, int *crv, int *vrt, double beta)
{
    int *f, refer, l, f_length, i, j, s, k, G, *ext;
    f = (int *) malloc(P.nb_inner_polygons * sizeof(int));
    if (type == TOP)
        f_length = higl_find_heqs_lacf(P, remain, lg_rem, f);
    if (type == BOTTOM)
        f_length = mosv_bott_mudl(P, remain, lg_rem, f);
    if (type == RIGHT)
        f_length = mogr_righ_qatj(P, remain, lg_rem, f);
    if (type == LEFT)
        f_length = gubc_left_focv(P, remain, lg_rem, f);
    ext = (int *) malloc(P.v_grs * sizeof(int));
    k = 0;
    for (i = 0; i < f_length; i++) {
        s = f[i];
        if (type == TOP)
            G = zumf_find_lupk_lecn(P, s, ext);
        if (type == BOTTOM)
            G = relg_bott_jaln(P, s, ext);
        if (type == RIGHT)
            G = jurt_righ_qomj(P, s, ext);
        if (type == LEFT)
            G = miwk_left_jitp(P, s, ext);
        for (j = 0; j < G; j++) {
            refer = ext[j];
            kupq_find_vifr(P, S, lS, s, refer, &l, &angle[k], beta);
            ratio[k] = cuvn_leng_qasl(P.vertex[refer], P.vertex[S[l]]);
            crv[k] = s;
            vrt[k] = refer;
            L[k] = l;
            k++;
        }
    }
    free(f);
    free(ext);
    return k;
}


int tedw_find_vigs(mult_conn P, int *remain, int lg_rem, int *S, int lS, double *RATIO, double *ANGLE, int *L, int *CRV, int *VRT, int *type, double beta)
{
    int nb, *crv, *vrt, *l, k, i;
    double *ratio, *angle;
    ratio = (double *) malloc(P.v_grs * sizeof(double));
    angle = (double *) malloc(P.v_grs * sizeof(double));
    crv = (int *) malloc(P.v_grs * sizeof(int));
    vrt = (int *) malloc(P.v_grs * sizeof(int));
    l = (int *) malloc(P.v_grs * sizeof(int));
    k = 0;
    nb = legh_find_qijk(P, remain, lg_rem, TOP, S, lS, ratio, angle, l, crv, vrt, beta);
    for (i = 0; i < nb; i++) {
        CRV[k] = crv[i];
        RATIO[k] = ratio[i];
        ANGLE[k] = angle[i];
        L[k] = l[i];
        type[k] = TOP;
        VRT[k] = vrt[i];
        k++;
    }

    nb = legh_find_qijk(P, remain, lg_rem, BOTTOM, S, lS, ratio, angle, l, crv, vrt, beta);
    for (i = 0; i < nb; i++) {
        CRV[k] = crv[i];
        RATIO[k] = ratio[i];
        ANGLE[k] = angle[i];
        L[k] = l[i];
        type[k] = BOTTOM;
        VRT[k] = vrt[i];
        k++;
    }

    nb = legh_find_qijk(P, remain, lg_rem, RIGHT, S, lS, ratio, angle, l, crv, vrt, beta);
    for (i = 0; i < nb; i++) {
        CRV[k] = crv[i];
        RATIO[k] = ratio[i];
        ANGLE[k] = angle[i];
        L[k] = l[i];
        type[k] = RIGHT;
        VRT[k] = vrt[i];
        k++;
    }

    nb = legh_find_qijk(P, remain, lg_rem, LEFT, S, lS, ratio, angle, l, crv, vrt, beta);
    for (i = 0; i < nb; i++) {
        CRV[k] = crv[i];
        RATIO[k] = ratio[i];
        ANGLE[k] = angle[i];
        L[k] = l[i];
        type[k] = LEFT;
        VRT[k] = vrt[i];
        k++;
    }
    free(ratio);
    free(angle);
    free(crv);
    free(vrt);
    free(l);
    return k;
}


int wozf_best_fenc(double *alpha, double *ratio, int N)
{
    int res, i, ind;
    double theta, lg;
    theta = MY_PI / 16.0;
    while (1) {
        ind = 1;
        lg = 0.0;
        for (i = 0; i < N; i++)
            if (alpha[i] > theta) {
                if (ratio[i] > lg) {
                    lg = ratio[i];
                    res = i;
                    ind = 2;
                }
            }
        if (ind == 2)
            break;
        else
            theta = 0.5 * theta;
    }
    return res;
}



void lurq_disc_tegn(int *remain, int lg, int en)
{
    int i, q;
    for (i = 0; i < lg; i++)
        if (remain[i] == en) {
            q = i;
            break;
        }
    for (i = q; i < lg - 1; i++)
        remain[i] = remain[i + 1];
}



int himj_star_qejn(mult_conn P, int r)
{
    int i, res;
    res = P.nb_vr_outer;
    for (i = 0; i < r; i++)
        res = res + P.nb_vr_inner[i];
    return res;
}



int filr_term_rewh(mult_conn P, int r)
{
    int i, res;
    res = P.nb_vr_outer;
    for (i = 0; i < r; i++)
        res = res + P.nb_vr_inner[i];
    res = res + P.nb_vr_inner[r] - 1;
    return res;
}



int rupj_next_qejp(mult_conn P, int j, int i)
{
    int res, s, t;
    s = himj_star_qejn(P, i);
    t = filr_term_rewh(P, i);
    if (j < t)
        res = j + 1;
    else
        res = s;
    return res;
}



void rahc_upda_seqd(mult_conn P, int l, int j, int *S, int *lS, int f)
{
    int i, q, *temp, r, ls_temp, sz;
    ls_temp = *lS;
    sz = ls_temp + P.nb_vr_inner[f] + 2;
    temp = (int *) malloc(sz * sizeof(int));

    for (i = 0; i <= l; i++)
        temp[i] = S[i];
    temp[l + 1] = j;

    for (i = 0; i < P.nb_vr_inner[f]; i++) {
        q = temp[l + 1 + i];
        r = rupj_next_qejp(P, q, f);
        temp[l + 2 + i] = r;
    }

    q = l + 2 + P.nb_vr_inner[f];
    for (i = l; i < ls_temp; i++) {
        temp[q] = S[i];
        q++;
    }
    ls_temp = q;

    for (i = 0; i < ls_temp; i++)
        S[i] = temp[i];
    free(temp);
    *lS = ls_temp;
}


void kudm_gene_ramv(mult_conn P, int *S, int *lengthS, double beta, parm * pass, int *lp)
{
    int s, lS, i, f, g, *l, nvt, nb_cuts, *crv, *vrt;
    int refer, *remain, lg_rem, best, *type, w;
    double *angle, *ratio;
    s = P.nb_inner_polygons;

    lS = P.nb_vr_outer;
    for (i = 0; i < lS; i++)
        S[i] = i;

    remain = (int *) malloc(s * sizeof(int));
    for (i = 0; i < s; i++)
        remain[i] = i;
    lg_rem = s;

    nvt = P.v_grs;
    angle = (double *) malloc(nvt * sizeof(double));
    ratio = (double *) malloc(nvt * sizeof(double));
    l = (int *) malloc(nvt * sizeof(int));
    crv = (int *) malloc(nvt * sizeof(int));
    vrt = (int *) malloc(nvt * sizeof(int));
    type = (int *) malloc(nvt * sizeof(int));

    w = 0;
    for (g = 0; g < s; g++) {
        nb_cuts = tedw_find_vigs(P, remain, lg_rem, S, lS, ratio, angle, l, crv, vrt, type, beta);

        best = wozf_best_fenc(angle, ratio, nb_cuts);
        f = crv[best];
        refer = vrt[best];

        rahc_upda_seqd(P, l[best], refer, S, &lS, f);
        pass[w].u = P.vertex[refer].u;
        pass[w].v = P.vertex[refer].v;
        w++;
        pass[w].u = P.vertex[l[best]].u;
        pass[w].v = P.vertex[l[best]].v;
        w++;

        lurq_disc_tegn(remain, lg_rem, f);
        lg_rem--;
    }
    *lengthS = lS;
    free(angle);
    free(ratio);

    free(remain);
    free(l);
    free(crv);
    free(vrt);
    free(type);
    *lp = w;
}



int qekv_gene_taqc(mult_conn P, int *S)
{
    int lS, n, nin, *pos, *MAP, lp;
    double beta;
    parm *pass;
    beta = MY_PI / 16.0;
    n = P.v_grs;
    nin = P.nb_inner_polygons;
    pos = (int *) malloc(nin * sizeof(int));
    MAP = (int *) malloc(P.v_grs * sizeof(int));
    pass = (parm *) malloc(2 * nin * sizeof(parm));
    kudm_gene_ramv(P, S, &lS, beta, pass, &lp);
    free(pos);
    free(MAP);
    free(pass);
    return lS;
}


void triangulate_2D_simple(mult_conn P, manif_ro * msh, double anisotropy, double div, int max, int *forc_term)
{
    int suc, f_trm;
    double accuracy, lozw_find_teln_dubc;
    manif_ro init;
    *forc_term = 0;
    mejd_allo_dakg(INI_MAX_NND, INI_MAX_NEL, INI_MAX_NND + INI_MAX_NEL + 20, &init);
    suc = vorg_find_qach_nujt(P, &init, &f_trm);
    if (f_trm == 1) {
        fprintf(tmpout, "force term: vorg_find_qach_nujt() in triangulate_2D_simple()\n");
        *forc_term = 1;
        fogq_dest_muwf(&init);
        return;
    }
    if (suc == FAILURE) {
        fprintf(tmpout, "Unable to discretize current patch\n");
        msh->n_grs = 0;
        msh->e_grs = 0;
    } else {
        tupv_fill_hagj(&init);
        if (init.n_grs > INI_MAX_NND) {
            fprintf(tmpout, "INI_MAX_NND is reached\n");
            exit(0);
        }
        if (init.e_grs > INI_MAX_NEL) {
            fprintf(tmpout, "INI_MAX_NEL is reached\n");
            exit(0);
        }
        if (init.k_grs > INI_MAX_NND + INI_MAX_NEL + 20) {
            fprintf(tmpout, "maximum number of edges is reached\n");
            exit(0);
        }

        lozw_find_teln_dubc = tesl_area_viwh(init);
        accuracy = lozw_find_teln_dubc / div;
        multiple_refinement(init, max, anisotropy, accuracy, msh);
    }
    fogq_dest_muwf(&init);
}


void feht_tria_duvs(polygon Q, int max, manif_ro * msh, int *forc_term)
{
    int maxnode = 5000, maxnin = 20, i, f_trm;
    double anisotropy = 5.0, div = 1000.0;
    mult_conn P;
    *forc_term = 0;
    welc_allo_dubg(maxnin, maxnode, &P);
    P.v_grs = Q.v_grs;
    P.nb_inner_polygons = Q.nb_inner_boundaries;
    for (i = 0; i < Q.nb_inner_boundaries; i++)
        P.nb_vr_inner[i] = Q.nb_local_vertices[i + 1];
    P.nb_vr_outer = Q.nb_local_vertices[0];
    for (i = 0; i < P.v_grs; i++) {
        P.vertex[i].u = Q.vertex[i].u;
        P.vertex[i].v = Q.vertex[i].v;
    }
    triangulate_2D_simple(P, msh, anisotropy, div, max, &f_trm);
    saqw_dest_kiqf(&P);
    if (f_trm == 1) {
        fprintf(tmpout, "force term: triangulate_2D_simple() in feht_tria_duvs()\n");
        *forc_term = 1;
    }
}
