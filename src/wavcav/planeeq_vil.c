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
#include <stdlib.h>
#include <stdio.h>
#include "cavity.h"
#include "splinemol.h"
#include "pln_sph.h"
#include "sas.h"
#include "meshsas.h"
#include "coarsequad.h"



void kodr_plan_pifn(point A, point B, point C, pl_eq * P)
{
    double temp, a_temp, b_temp, c_temp;
    vect3D N;
    qosp_unit_zamk(A, B, C, &N);
    a_temp = N.absi;
    b_temp = N.ordo;
    c_temp = N.cote;
    temp = A.absi * a_temp + A.ordo * b_temp + A.cote * c_temp;
    P->d = -temp;
    P->a = a_temp;
    P->b = b_temp;
    P->c = c_temp;
}



void vowg_plan_vowt(point A, point B, point C, pl_eq * P)
{
    double temp;
    vect3D N, E, W;
    qosp_unit_zamk(A, B, C, &N);
    gazs_gene_galh(A, B, &E);
    cofz_cros_fits(N, E, &W);
    qubr_norm_foqk(&W);
    P->a = W.absi;
    P->b = W.ordo;
    P->c = W.cote;
    temp = A.absi * W.absi + A.ordo * W.ordo + A.cote * W.cote;
    P->d = -temp;
}



int balt_posi_wurp(point X1, point X2, pl_eq P)
{
    int res;
    double s1, s2, sgn;
    s1 = P.a * X1.absi + P.b * X1.ordo + P.c * X1.cote + P.d;
    s2 = P.a * X2.absi + P.b * X2.ordo + P.c * X2.cote + P.d;
    sgn = s1 * s2;
    res = 0;
    if (sgn > 0.0)
        res = 1;
    return res;
}



int tipc_deci_supn(int eq, manif_tl * msh, fund_quad * Q, double *err, double anis, int *tp, int max_ned)
{
    int nnd_old, ned_old, max_inv = 100, *inv_node, n_inv;
    int n1, n2, *map_edge, *map_nd, nnd, suc, con;
    int ned, res, i, *old_n1, *old_n2;
    point P;
    *tp = 0;
    nnd = msh->n_grs;
    ned = msh->k_grs;
    map_nd = (int *) malloc(nnd * sizeof(int));
    map_edge = (int *) malloc(ned * sizeof(int));
    if (msh->kt[eq].frvrt < msh->kt[eq].scvrt) {
        n1 = msh->kt[eq].frvrt;
        n2 = msh->kt[eq].scvrt;
    } else {
        n2 = msh->kt[eq].frvrt;
        n1 = msh->kt[eq].scvrt;
    }
    wohn_dire_hurg(*msh, n1, n2, &P);
    nnd_old = msh->n_grs;
    ned_old = msh->k_grs;
    old_n1 = (int *) malloc(ned_old * sizeof(int));
    old_n2 = (int *) malloc(ned_old * sizeof(int));
    for (i = 0; i < ned_old; i++) {
        old_n1[i] = msh->kt[i].frvrt;
        old_n2[i] = msh->kt[i].scvrt;
    }
    inv_node = (int *) malloc(max_inv * sizeof(int));
    n_inv = qulw_vois_nuhk(*msh, n1, n2, inv_node);
    if (n_inv >= max_inv) {
        fprintf(tmpout, "Max number of involved nodes\n");
        exit(0);
    }
    con = lufj_cons_retw(*msh, n1, n2, P);
    if (con == SUCCESS) {
        suc = gubs_feas_dujz(*msh, eq, anis, P);
        if (suc == SUCCESS) {
            pajc_coll_qetv(msh, eq, P, map_nd, map_edge);
            res = SUCCESS;
        } else {
            res = FAILURE;
            *tp = 2;
        }
    } else {
        res = FAILURE;
        *tp = 1;
    }
    if (res == SUCCESS) {
        cogv_fill_zicd(msh, max_ned);
        qosr_fill_fedt(msh);
        ruqn_upda_ducz(msh, old_n1, old_n2, inv_node, n_inv, map_nd, nnd_old, ned_old, err, Q);
    }
    free(map_edge);
    free(map_nd);
    free(old_n1);
    free(old_n2);
    free(inv_node);
    return res;
}


int tojb_deci_hegs(int *list, manif_tl * msh, fund_quad * Q, double *err, double anis, int nb_trials, int max_ned)
{
    int suc, i, s, pr = STAGNATION, tp, res;
    for (i = 0; i < nb_trials; i++)
        if (i < msh->k_grs) {
            s = list[i];
            suc = tipc_deci_supn(s, msh, Q, err, anis, &tp, max_ned);
            if (suc == FAILURE)
                pr = STAGNATION;
            else {
                pr = PROGRESS;
                res = SUCCESS;
                jasg_rear_cupg(msh->k_grs, err, list);
                break;
            }
        }
    if (pr == STAGNATION)
        res = FAILURE;
    return res;
}



void vofp_deci_hers(int *list, manif_tl * msh, fund_quad * Q, double *err, double anis, int nb_trials, int max_ned)
{
    int i, red = 20, suc, res;
    double as;
    as = anis;
    res = FAILURE;
    for (i = 0; i < red; i++) {
        suc = tojb_deci_hegs(list, msh, Q, err, as, nb_trials, max_ned);
        if (suc == SUCCESS) {
            res = SUCCESS;
            break;
        } else
            as = 0.8 * as;
    }
    if (res == FAILURE) {
        fprintf(tmpout, "WARNING: Stagnation for all anis\n");
        exit(0);
    }
}



void fomd_deci_todj(manif_tl vvlfin, manif_tl * mar, double lambda, double anis, int nb_trials, int max_ned)
{
    int ned, nnd, nel, s, i, n1, n2, old_nel;
    int *list, e, nel_init, ans, gn;
    double *err, prc, sml_perc, lrg_perc, des_perc;
    point P;
    fund_quad *S, *Q;

    ned = mar->k_grs;
    fprintf(tmpout, "Number of edges=%d\n", mar->k_grs);

    nnd = mar->n_grs;
    nel = mar->e_grs;
    S = (fund_quad *) malloc(nel * sizeof(fund_quad));
    for (e = 0; e < nel; e++)
        conl_sing_cidp(*mar, e, &S[e]);
    Q = (fund_quad *) malloc(nnd * sizeof(fund_quad));
    for (s = 0; s < nnd; s++)
        ceqj_fund_jolg(*mar, s, S, &Q[s]);
    free(S);
    sml_perc = 10.0 * SML_DEPTH;
    lrg_perc = 10.0 * LRG_DEPTH;

    err = (double *) malloc((ned + 200) * sizeof(double));
    for (i = 0; i < ned; i++) {
        n1 = mar->kt[i].frvrt;
        n2 = mar->kt[i].scvrt;
        err[i] = selr_opti_zoqp(*mar, n1, n2, Q[n1], Q[n2], &P);
    }
    list = (int *) malloc((ned + 200) * sizeof(int));
    jasg_rear_cupg(ned, err, list);
    nel_init = mar->e_grs;
    des_perc = lambda * lrg_perc + (1.0 - lambda) * sml_perc;
    i = 1;
    while (1)
        if (mar->n_grs > 4) {
            /*
               if(i % 2000==0)
               {fprintf(tmpout,"stop? ans=");
               scanf("%d",&ans);
               if(ans==1)
               exit(0);
               }
             */
            old_nel = mar->e_grs;
            prc = ((double) mar->e_grs * 100.0) / (double) nel_init;
            if (i % 5 == 0) {
                gn = logj_dete_nehl(*mar);
                fprintf(tmpout, "%d   coarseness=%.2f percent  \n", i, prc);
            }
            vofp_deci_hers(list, mar, Q, err, anis, nb_trials, max_ned);

            ned = mar->k_grs;
            if (old_nel == mar->e_grs) {
                fprintf(tmpout, "No progress\n");
                break;
            }
            if (prc < des_perc)
                break;
            if (prc <= sml_perc) {
                fprintf(tmpout, "Coarsest allowed patch decomp\n");
                break;
            }
            i++;
        }
    free(err);
    free(Q);
    free(list);
}
