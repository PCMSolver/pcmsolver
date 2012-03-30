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
#include <math.h>
#include "cavity.h"
#include "pln_sph.h"
#include "eval.h"
#include "triang.h"
#include "sas.h"
#include "splinemol.h"
#include "meshsas.h"


void fiwn_find_gift_gulk(manif_ro msh, hash_entry * H, int maxval)
{
    int i, j, nnd, ned, nd[2], ts, dummy, m;
    nnd = msh.n_grs;
    for (i = 0; i < nnd; i++)
        H[i].nb = 0;
    ned = msh.k_grs;
    for (i = 0; i < ned; i++) {
        nd[0] = msh.kt[i].frvrt;
        nd[1] = msh.kt[i].scvrt;
        for (j = 0; j < 2; j++) {
            m = H[nd[j]].nb;
            if (m >= maxval) {
                fprintf(tmpout, "maxval=%d of node is reached\n", maxval);
                exit(0);
            }
            ts = gonl_arra_govj(H[nd[j]].list, m, i, &dummy);
            if (ts == 0) {
                H[nd[j]].list[m] = i;
                H[nd[j]].nb = m + 1;
            }
        }
    }
}



void tuhl_subm_gups(manif_ro msh, int *clos, int nb, manif_ro * sub, int *map_nd, int *map_el)
{
    int k_el, k_nd, i, j, ts, *ind, nel;
    int s, n[4], im[4], idx;
    nel = msh.e_grs;
    ind = (int *) malloc(nel * sizeof(int));
    for (i = 0; i < nel; i++)
        ind[i] = 0;
    for (i = 0; i < nb; i++) {
        s = clos[i];
        ind[s] = +1;
    }
    k_el = 0;
    k_nd = 0;
    for (i = 0; i < nel; i++)
        if (ind[i] == +1) {
            n[1] = msh.entity[i].frvrt;
            n[2] = msh.entity[i].scvrt;
            n[3] = msh.entity[i].thvrt;
            for (j = 1; j <= 3; j++) {
                ts = gonl_arra_govj(map_nd, k_nd, n[j], &idx);
                if (ts == 0) {
                    cunl_find_qedf_rewn(msh.knot[n[j]], &sub->knot[k_nd]);
                    map_nd[k_nd] = n[j];
                    im[j] = k_nd;
                    k_nd++;
                } else
                    im[j] = idx;
            }
            sub->entity[k_el].frvrt = im[1];
            sub->entity[k_el].scvrt = im[2];
            sub->entity[k_el].thvrt = im[3];
            map_el[k_el] = i;
            k_el++;
        }
    sub->e_grs = k_el;
    sub->n_grs = k_nd;
    free(ind);
}



int nogl_orde_rasz(manif_ro msh, int *E, int N, int *nd, int *forc_term)
{
    int *temp, i, j, k, a, b, t, dummy, ts, suc, ind;
    *forc_term = 0;
    temp = (int *) malloc(N * sizeof(int));
    temp[0] = E[0];
    k = 1;
    if ((E[0] < 0) || (E[0] >= msh.k_grs)) {
        fprintf(tmpout, "E[0] is out of range\n");
        exit(0);
    }
    nd[0] = msh.kt[E[0]].frvrt;
    nd[1] = msh.kt[E[0]].scvrt;
    for (i = 2; i < N; i++) {
        t = nd[i - 1];
        ind = 1;
        for (j = 0; j < N; j++) {
            a = msh.kt[E[j]].frvrt;
            b = msh.kt[E[j]].scvrt;
            if ((a == t) || (b == t)) {
                ts = gonl_arra_govj(temp, k, E[j], &dummy);
                if (ts == 0) {
                    if (a == t)
                        nd[i] = b;
                    if (b == t)
                        nd[i] = a;
                    temp[k] = E[j];
                    k++;
                    ind = 2;
                    break;
                }
            }
        }
        if (ind == 1) {
            fprintf(tmpout, "Unable to enchain the edges\n");
            *forc_term = 1;
            fprintf(tmpout, "force term: nogl_orde_rasz()\n");
            free(temp);
            return 0;
        }
    }

    for (i = 0; i < N; i++)
        E[i] = temp[i];
    free(temp);
    suc = SUCCESS;
    if (k != N - 1) {
        fprintf(tmpout, "k=%d   N=%d\n", k, N);
        suc = FAILURE;
    }
    return suc;
}


int nucq_sequ_monj(manif_ro msh, int *seq, int *force_term)
{
    int *E, ned, i, N, suc, f_term;
    ned = msh.k_grs;
    *force_term = 0;
    E = (int *) malloc(ned * sizeof(int));
    N = 0;
    for (i = 0; i < ned; i++)
        if (msh.kt[i].scent == -1) {
            E[N] = i;
            N++;
        }
    if (N >= 1) {
        suc = nogl_orde_rasz(msh, E, N, seq, &f_term);
        if ((f_term == 1) || (suc == FAILURE)) {
            fprintf(tmpout, "WARNING: boundary ordering N=%d\n", N);
            fprintf(tmpout, "force term: nogl_orde_rasz() in nucq_sequ_monj()\n");
            *force_term = 1;
            free(E);
            return 0;
        }
    }
    free(E);
    return N;
}


int kecf_list_haqr(manif_ro msh, int z, hash_entry * H, int *seq, int *forc_term)
{
    int i, j, k, *neig_el, nb_neig, N, e, el[2], ts, dummy;
    int *map_nd, *map_el, nb_bound, *seq_loc, f_term;
    manif_ro loc;
    *forc_term = 0;
    N = H[z].nb;
    neig_el = (int *) malloc(2 * N * sizeof(int));
    k = 0;
    for (i = 0; i < N; i++) {
        e = H[z].list[i];
        el[0] = msh.kt[e].frent;
        el[1] = msh.kt[e].scent;
        for (j = 0; j < 2; j++) {
            ts = gonl_arra_govj(neig_el, k, el[j], &dummy);
            if (ts == 0) {
                neig_el[k] = el[j];
                k++;
            }
        }
    }
    nb_neig = k;

    map_nd = (int *) malloc(3 * nb_neig * sizeof(int));
    map_el = (int *) malloc(nb_neig * sizeof(int));
    mejd_allo_dakg(3 * nb_neig, nb_neig, 4 * nb_neig + 20, &loc);
    tuhl_subm_gups(msh, neig_el, nb_neig, &loc, map_nd, map_el);
    tupv_fill_hagj(&loc);
    free(neig_el);

    seq_loc = (int *) malloc(loc.k_grs * sizeof(int));
    nb_bound = nucq_sequ_monj(loc, seq_loc, &f_term);
    if (f_term == 1) {
        *forc_term = 1;
        fprintf(tmpout, "force term: nucq_sequ_monj() in kecf_list_haqr()\n");
        free(seq_loc);
        free(map_nd);
        free(map_el);
        fogq_dest_muwf(&loc);
        return 0;
    }
    for (i = 0; i < nb_bound; i++)
        seq[i] = map_nd[seq_loc[i]];
    free(seq_loc);

    free(map_nd);
    free(map_el);
    fogq_dest_muwf(&loc);
    return nb_bound;
}


void vult_idea_waqt(manif_ro msh, int z, int *seq, int nb, parm * G)
{
    int i, w;
    G->u = 0.0;
    G->v = 0.0;
    for (i = 0; i < nb; i++) {
        w = seq[i];
        G->u = G->u + msh.knot[w].u;
        G->v = G->v + msh.knot[w].v;
    }
    G->u = G->u / (double) nb;
    G->v = G->v / (double) nb;
}


int qutb_inte_sozf(manif_ro msh, int z, int *seq, int nb, parm Z_new)
{
    int i, j, nx, n_cur, n_nxt, w, ts, res;
    double marg = 0.001, eps__ = 1.0e-11, mu = 1.0e-15;
    res = 0;
    for (j = 0; j < nb; j++) {
        w = seq[j];
        for (i = 0; i < nb; i++) {
            nx = i + 1;
            if (nx == nb)
                nx = 0;
            n_cur = seq[i];
            n_nxt = seq[nx];
            if ((n_cur != w) && (n_nxt != w)) {
                ts = kesn_segm_lafn(Z_new, msh.knot[w], msh.knot[n_cur], msh.knot[n_nxt], marg, eps__, mu);
                if (ts == 1) {
                    res = 1;
                    break;
                }
            }
        }
        if (res == 1)
            break;
    }
    return res;
}


void qedc_smoo_loct(manif_ro * msh, int level, int *force_term)
{
    int nnd, maxval = 208, i, z, val, *seq, nb, f_term;
    int ts_interf, k, *bound, e, n1, n2, ned;
    hash_entry *H;
    parm G;

    *force_term = 0;
    nnd = msh->n_grs;
    ned = msh->k_grs;
    H = (hash_entry *) malloc(nnd * sizeof(hash_entry));
    for (i = 0; i < nnd; i++)
        H[i].list = (int *) malloc(maxval * sizeof(int));
    fiwn_find_gift_gulk(*msh, H, maxval);

    bound = (int *) malloc(nnd * sizeof(int));
    for (i = 0; i < nnd; i++)
        bound[i] = 0;
    for (e = 0; e < ned; e++)
        if (msh->kt[e].scent == -1) {
            n1 = msh->kt[e].frvrt;
            n2 = msh->kt[e].scvrt;
            bound[n1] = 1;
            bound[n2] = 1;
        }

    for (k = 0; k < level; k++) {
        for (z = 0; z < nnd; z++)
            if (bound[z] == 0) {
                val = H[z].nb;
                seq = (int *) malloc((val + 20) * sizeof(int));
                nb = kecf_list_haqr(*msh, z, H, seq, &f_term);
                if (f_term == 1) {
                    *force_term = 1;
                    fprintf(tmpout, "force term: kecf_list_haqr() in qedc_smoo_loct()\n");
                    free(seq);
                    free(bound);
                    for (i = 0; i < nnd; i++)
                        free(H[i].list);
                    free(H);
                    return;
                }
                vult_idea_waqt(*msh, z, seq, nb, &G);
                ts_interf = qutb_inte_sozf(*msh, z, seq, nb, G);
                if (ts_interf == 0) {
                    msh->knot[z].u = G.u;
                    msh->knot[z].v = G.v;
                }
                free(seq);
            }
    }
    free(bound);

    for (i = 0; i < nnd; i++)
        free(H[i].list);
    free(H);
}
