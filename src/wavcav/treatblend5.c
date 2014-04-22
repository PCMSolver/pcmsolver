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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cavity.h"
#include "pln_sph.h"
#include "eval.h"
#include "triang.h"
#include "sas.h"
#include "partcas.h"
#include "splinemol.h"
#include "coarsequad.h"
#include "geodesic.h"
#include "meshsas.h"


int kuts_next_kolv(telolf T, int s)
{
    int res = -1;
    if (T.frvrt == s)
        res = T.scvrt;
    else if (T.scvrt == s)
        res = T.thvrt;
    else if (T.thvrt == s)
        res = T.frvrt;
    return res;
}



int kotc_test_kolj(manif_tl msh, int p, int q)
{
    int *nd_p, *nd_q, *cm, k, i, res;
    int ts, dummy, nx_p, nx_q;
    nd_p = (int *) malloc(3 * sizeof(int));
    nd_q = (int *) malloc(3 * sizeof(int));
    cm = (int *) malloc(3 * sizeof(int));
    nd_p[0] = msh.entity[p].frvrt;
    nd_p[1] = msh.entity[p].scvrt;
    nd_p[2] = msh.entity[p].thvrt;

    nd_q[0] = msh.entity[q].frvrt;
    nd_q[1] = msh.entity[q].scvrt;
    nd_q[2] = msh.entity[q].thvrt;
    k = 0;
    for (i = 0; i < 3; i++) {
        ts = gonl_arra_govj(nd_q, 3, nd_p[i], &dummy);
        if (ts == 1) {
            cm[k] = nd_p[i];
            k++;
        }
    }
    if (k != 2) {
        fprintf(tmpout, "Elements must share two nodes\n");
        fprintf(tmpout, "[p,q]=[%d,%d]\n", p, q);
        exit(0);
    }
    res = 1;
    for (i = 0; i < 2; i++) {
        nx_p = kuts_next_kolv(msh.entity[p], cm[i]);
        nx_q = kuts_next_kolv(msh.entity[q], cm[i]);
        if ((nx_p == -1) || (nx_q == -1)) {
            fprintf(tmpout, "Unable to find next node\n");
            exit(0);
        }
        if (nx_p == nx_q)
            res = 0;
    }
    free(nd_p);
    free(nd_q);
    free(cm);
    return res;
}



int fevd_loca_logh(manif_tl msh, int s)
{
    int ts = 1, *nd, ed[3], i, e, n1, n2, tr1, tr2, dummy;
    nd = (int *) malloc(3 * sizeof(int));
    nd[0] = msh.entity[s].frvrt;
    nd[1] = msh.entity[s].scvrt;
    nd[2] = msh.entity[s].thvrt;
    ed[0] = msh.entity[s].frkt;
    ed[1] = msh.entity[s].sckt;
    ed[2] = msh.entity[s].trkt;
    if ((ed[0] < 0) || (ed[1] < 0) || (ed[2] < 0)) {
        fprintf(tmpout, "edges of el[%d]=[%d,%d,%d]\n", s, ed[0], ed[1], ed[2]);
        exit(0);
    }

    for (i = 0; i < 3; i++) {
        e = ed[i];
        n1 = msh.kt[e].frvrt;
        n2 = msh.kt[e].scvrt;
        tr1 = gonl_arra_govj(nd, 3, n1, &dummy);
        tr2 = gonl_arra_govj(nd, 3, n2, &dummy);
        if ((tr1 == 0) || (tr2 == 0)) {
            ts = 0;
            break;
        }
    }
    free(nd);
    if (ts == 0) {
        fprintf(tmpout, "element=%d with bad integrity1\n", s);
        exit(0);
    }

    for (i = 0; i < 3; i++) {
        e = ed[i];
        if ((msh.kt[e].frent != s) && (msh.kt[e].scent != s)) {
            ts = 0;
            break;
        }
    }
    if (ts == 0) {
        fprintf(tmpout, "nodes=[%d,%d,%d]\n", msh.entity[s].frvrt, msh.entity[s].scvrt, msh.entity[s].thvrt);
        fprintf(tmpout, "edges=[%d,%d,%d]\n", msh.entity[s].frkt, msh.entity[s].sckt, msh.entity[s].trkt);
        fprintf(tmpout, "element=%d with bad integrity2\n", s);
        exit(0);
    }
    return ts;
}



int local_element_integrity_forc(manif_tl msh, int s)
{
    int ts = 1, *nd, ed[3], i, e, n1, n2, tr1, tr2, dummy;
    nd = (int *) malloc(3 * sizeof(int));
    nd[0] = msh.entity[s].frvrt;
    nd[1] = msh.entity[s].scvrt;
    nd[2] = msh.entity[s].thvrt;
    ed[0] = msh.entity[s].frkt;
    ed[1] = msh.entity[s].sckt;
    ed[2] = msh.entity[s].trkt;
    if ((ed[0] < 0) || (ed[1] < 0) || (ed[2] < 0)) {
        fprintf(tmpout, "edges of el[%d]=[%d,%d,%d]\n", s, ed[0], ed[1], ed[2]);
        free(nd);
        return 0;
    }

    for (i = 0; i < 3; i++) {
        e = ed[i];
        n1 = msh.kt[e].frvrt;
        n2 = msh.kt[e].scvrt;
        tr1 = gonl_arra_govj(nd, 3, n1, &dummy);
        tr2 = gonl_arra_govj(nd, 3, n2, &dummy);
        if ((tr1 == 0) || (tr2 == 0)) {
            ts = 0;
            break;
        }
    }
    free(nd);
    if (ts == 0) {
        fprintf(tmpout, "element=%d with bad integrity1\n", s);
        return 0;
    }

    for (i = 0; i < 3; i++) {
        e = ed[i];
        if ((msh.kt[e].frent != s) && (msh.kt[e].scent != s)) {
            ts = 0;
            break;
        }
    }
    if (ts == 0) {
        fprintf(tmpout, "nodes=[%d,%d,%d]\n", msh.entity[s].frvrt, msh.entity[s].scvrt, msh.entity[s].thvrt);
        fprintf(tmpout, "edges=[%d,%d,%d]\n", msh.entity[s].frkt, msh.entity[s].sckt, msh.entity[s].trkt);
        fprintf(tmpout, "element=%d with bad integrity2\n", s);
        return 0;
    }
    return 1;
}



int qult_loca_qulm(manif_tl msh, int s)
{
    int ts = 1, e[2], i, f;
    e[0] = msh.kt[s].frent;
    e[1] = msh.kt[s].scent;
    for (i = 0; i < 2; i++) {
        f = e[i];
        if ((msh.entity[f].frkt != s) && (msh.entity[f].sckt != s) && (msh.entity[f].trkt != s)) {
            ts = 0;
            break;
        }
    }
    if (ts == 0) {
        fprintf(tmpout, "kt=%d with bad integrity1\n", s);
        fprintf(tmpout, "element1=%d  incident edges=[%d,%d,%d]\n", e[0], msh.entity[e[0]].frkt, msh.entity[e[0]].sckt, msh.entity[e[0]].trkt);
        fprintf(tmpout, "element2=%d  incident edges=[%d,%d,%d]\n", e[1], msh.entity[e[1]].frkt, msh.entity[e[1]].sckt, msh.entity[e[1]].trkt);
        exit(0);
    }
    return ts;
}



int local_edge_integrity_forc(manif_tl msh, int s)
{
    int ts = 1, e[2], i, f;
    e[0] = msh.kt[s].frent;
    e[1] = msh.kt[s].scent;
    for (i = 0; i < 2; i++) {
        f = e[i];
        if ((msh.entity[f].frkt != s) && (msh.entity[f].sckt != s) && (msh.entity[f].trkt != s)) {
            ts = 0;
            break;
        }
    }
    if (ts == 0) {
        fprintf(tmpout, "kt=%d with bad integrity1\n", s);
        fprintf(tmpout, "element1=%d  incident edges=[%d,%d,%d]\n", e[0], msh.entity[e[0]].frkt, msh.entity[e[0]].sckt, msh.entity[e[0]].trkt);
        fprintf(tmpout, "element2=%d  incident edges=[%d,%d,%d]\n", e[1], msh.entity[e[1]].frkt, msh.entity[e[1]].sckt, msh.entity[e[1]].trkt);
        return 0;
    }
    return 1;
}


void hajk_chec_wizn(manif_tl msh)
{
    int nnd, nel, ned, s;

    nnd = msh.n_grs;
    nel = msh.e_grs;
    ned = msh.k_grs;

    for (s = 0; s < nel; s++)
        fevd_loca_logh(msh, s);

    for (s = 0; s < ned; s++)
        qult_loca_qulm(msh, s);
    fprintf(tmpout, "GOOD MANIFOLD INTEGRITY\n");
}


int check_mesh_integrity_forc(manif_tl msh)
{
    int nnd, nel, ned, s, res, ts;

    nnd = msh.n_grs;
    nel = msh.e_grs;
    ned = msh.k_grs;

    res = 1;
    for (s = 0; s < nel; s++) {
        ts = local_element_integrity_forc(msh, s);
        if (ts == 0) {
            res = 0;
            break;
        }
    }

    for (s = 0; s < ned; s++) {
        ts = local_edge_integrity_forc(msh, s);
        if (ts == 0) {
            res = 0;
            break;
        }
    }
    if (res == 1)
        fprintf(tmpout, "GOOD MANIFOLD INTEGRITY\n");
    return res;
}



int qunf_expa_foqg(manif_tl * msh, int *bulk, hash_entry * H)
{
    int suc = FAILURE, el, nel, id, ort, n2, n3, nb;
    nel = msh->e_grs;
    nb = 0;
    for (el = 0; el < nel; el++)
        if (bulk[el] == -1) {
            id = vist_find_jetl(*msh, el, bulk, H);
            if (id != -1) {
                ort = kotc_test_kolj(*msh, el, id);
                if (ort == 0) {
                    n2 = msh->entity[el].scvrt;
                    n3 = msh->entity[el].thvrt;
                    msh->entity[el].scvrt = n3;
                    msh->entity[el].thvrt = n2;
                }
                bulk[el] = +1;
                nb++;
                suc = SUCCESS;
            }
        }
    return suc;
}


void kumn_inci_fuzr(manif_tl msh, hash_entry * H)
{
    int nel, i, j, k, e[3], E1, E2;
    nel = msh.e_grs;
    for (i = 0; i < nel; i++) {
        e[0] = msh.entity[i].frkt;
        e[1] = msh.entity[i].sckt;
        e[2] = msh.entity[i].trkt;
        k = 0;
        for (j = 0; j < 3; j++) {
            E1 = msh.kt[e[j]].frent;
            E2 = msh.kt[e[j]].scent;
            if (E1 == i) {
                H[i].list[k] = E2;
                k++;
            } else if (E2 == i) {
                H[i].list[k] = E1;
                k++;
            }
        }
        H[i].nb = k;
    }
    for (i = 0; i < nel; i++)
        if (H[i].nb != 3)
            fprintf(tmpout, "WARNING: manifold might not be closed\n");
}



int rukn_mani_jotq(manif_tl * msh, hash_entry * H)
{
    int *bulk, i, suc, seed, q, nb_comp_max = 4;
    int nb_comp, ind, nel, ned, p;
    nel = msh->e_grs;
    bulk = (int *) malloc(nel * sizeof(int));
    nb_comp = 1;
    seed = 0;
    for (i = 0; i < nel; i++)
        if (i != seed)
            bulk[i] = -1;
    for (q = 0; q < nb_comp_max; q++) {
        bulk[seed] = +1;
        ned = msh->k_grs;
        for (p = 0; p < ned; p++) {
            suc = qunf_expa_foqg(msh, bulk, H);
            if (suc == FAILURE)
                break;
        }
        ind = 1;
        for (i = 0; i < nel; i++)
            if (bulk[i] == -1) {
                seed = i;
                ind = 2;
                break;
            }
        if ((ind == 2) && (q == nb_comp_max)) {
            fprintf(tmpout, "Max nb manifold components is reached\n");
            exit(0);
        }
        if (ind == 1)
            break;
        nb_comp++;
    }
    fprintf(tmpout, "Complete expansion\n");
    free(bulk);
    return nb_comp;
}



void kecq_orie_watk(manif_tl * msh, int max_ned)
{
    int nel, i, nb_comp, ort, e1, e2;
    hash_entry *H;
    nel = msh->e_grs;
    H = (hash_entry *) malloc(nel * sizeof(hash_entry));
    for (i = 0; i < nel; i++)
        H[i].list = (int *) malloc(3 * sizeof(int));
    kumn_inci_fuzr(*msh, H);
    nb_comp = rukn_mani_jotq(msh, H);
    fprintf(tmpout, "NB COMPONENTS=%d\n", nb_comp);
    cogv_fill_zicd(msh, max_ned);
    for (i = 0; i < msh->k_grs; i++) {
        e1 = msh->kt[i].frent;
        e2 = msh->kt[i].scent;
        if (e2 != -1) {
            ort = kotc_test_kolj(*msh, e1, e2);
            if (ort == 0) {
                fprintf(tmpout, "kt=%d\n", i);
                fprintf(tmpout, "WARNING: Inconsistent orientation\n");
                fprintf(tmpout, "element[%d]=[%d,%d,%d]   ", e1, msh->entity[e1].frvrt, msh->entity[e1].scvrt, msh->entity[e1].thvrt);
                fprintf(tmpout, "element[%d]=[%d,%d,%d]\n", e2, msh->entity[e2].frvrt, msh->entity[e2].scvrt, msh->entity[e2].thvrt);
                exit(0);
            }
        }
    }
    for (i = 0; i < nel; i++)
        free(H[i].list);
    free(H);
    fprintf(tmpout, "GOOD LOCAL ORIENTATIONS\n");
}


int sutb_find_zokw_codh(trmsrf surf)
{
    int nb, nin, i;
    nb = surf.cc.N;
    nin = surf.nb_inner;
    for (i = 0; i < nin; i++)
        nb = nb + surf.inner[i].N;
    return nb;
}


void jizw_allo_vipg(trmsrf * surf1, int nb_surf1, int nb_surf2, int nb_surf3, int mx_edge_m, prat_main_m * GM)
{
    int i, k, nb_msh, nb_c;
    nb_msh = nb_surf1 + nb_surf2 + nb_surf3;
    GM->N = (sp_nd_m *) malloc(nb_msh * sizeof(sp_nd_m));
    k = 0;
    for (i = 0; i < nb_surf2; i++) {
        GM->N[k].corner = (int *) malloc(4 * sizeof(int));
        GM->N[k].inc_edge_m = (int *) malloc(4 * sizeof(int));
        k++;
    }
    for (i = 0; i < nb_surf3; i++) {
        GM->N[k].corner = (int *) malloc(3 * sizeof(int));
        GM->N[k].inc_edge_m = (int *) malloc(3 * sizeof(int));
        k++;
    }
    for (i = 0; i < nb_surf1; i++) {
        nb_c = sutb_find_zokw_codh(surf1[i]);
        GM->N[k].corner = (int *) malloc(nb_c * sizeof(int));
        GM->N[k].inc_edge_m = (int *) malloc(nb_c * sizeof(int));
        k++;
    }
    GM->E = (kt_m *) malloc(mx_edge_m * sizeof(kt_m));
}


void kulf_dest_lajs(int nb_surf1, int nb_surf2, int nb_surf3, int mx_edge_m, prat_main_m * GM)
{
    int i, nb_msh;
    nb_msh = nb_surf1 + nb_surf2 + nb_surf3;
    for (i = 0; i < nb_msh; i++)
        free(GM->N[i].corner);
    free(GM->N);
    free(GM->E);
}


void najt_prop_pics(megamanif MG, int *NND, int *NEL, int *NED)
{
    int N, nnd, nel, ned, i;
    N = MG.mw_grs;
    nnd = 0;
    nel = 0;
    ned = 0;
    for (i = 0; i < N; i++) {
        nnd = nnd + MG.msh[i].n_grs;
        nel = nel + MG.msh[i].e_grs;
        ned = ned + MG.msh[i].k_grs;
    }
    *NND = nnd;
    *NEL = nel;
    *NED = ned;
}



void dulw_find_tigv_dubr(manif_tl msh)
{
    int nnd, i;
    FILE *fp;
    fp = fopen("nodes_c.dat", "w");
    nnd = msh.n_grs;
    for (i = 0; i < nnd; i++)
        fprintf(fp, "%f  %f   %f\n", msh.knot[i].absi, msh.knot[i].ordo, msh.knot[i].cote);
    fclose(fp);
}


void pows_popu_peql(manif_tl msh, manif_tl * vvlh)
{
    int nnd, i;
    nnd = msh.n_grs;
    for (i = 0; i < nnd; i++) {
        vvlh->knot[i].absi = msh.knot[i].absi;
        vvlh->knot[i].ordo = msh.knot[i].ordo;
        vvlh->knot[i].cote = msh.knot[i].cote;
    }
    vvlh->n_grs = nnd;
}


void deqt_find_vegc_cobn(manif_tl msh)
{
    int nel, i;
    FILE *fp;
    fp = fopen("elements_c.dat", "w");
    nel = msh.e_grs;
    for (i = 0; i < nel; i++)
        fprintf(fp, "%d  %d  %d\n", msh.entity[i].frvrt, msh.entity[i].scvrt, msh.entity[i].thvrt);
    fclose(fp);
}


void hisl_popu_nuvw(manif_tl msh, manif_tl * vvlh)
{
    int nel, i;
    nel = msh.e_grs;
    for (i = 0; i < nel; i++) {
        vvlh->entity[i].frvrt = msh.entity[i].frvrt;
        vvlh->entity[i].scvrt = msh.entity[i].scvrt;
        vvlh->entity[i].thvrt = msh.entity[i].thvrt;
    }
    vvlh->e_grs = nel;
}


void wegf_find_hord_luvk(manif_tl msh)
{
    int nnd, i;
    FILE *fp;
    fp = fopen("nodes_f.dat", "w");
    nnd = msh.n_grs;
    for (i = 0; i < nnd; i++)
        fprintf(fp, "%.20f  %.20f   %.20f\n", msh.knot[i].absi, msh.knot[i].ordo, msh.knot[i].cote);
    fclose(fp);
}


void somn_popu_wolm(manif_tl msh, manif_tl * msh_f)
{
    int nnd, i;
    nnd = msh.n_grs;
    for (i = 0; i < nnd; i++) {
        msh_f->knot[i].absi = msh.knot[i].absi;
        msh_f->knot[i].ordo = msh.knot[i].ordo;
        msh_f->knot[i].cote = msh.knot[i].cote;
    }
    msh_f->n_grs = msh.n_grs;
}


void somt_find_qoph_hegd(manif_tl msh, rgb_lk * col, trmsrf * surf1, trmsrf * surf2, trmsrf * surf3, supp_surf * ssr)
{
    int nel, i, tp, id, w, n1, nb;
    double x, y, z, r, dis, diff, err;
    FILE *fp;
    fp = fopen("elements_f.dat", "w");
    nel = msh.e_grs;
    err = 0.0;
    nb = 0;
    for (i = 0; i < nel; i++) {
        tp = ssr[i].s_type;
        id = ssr[i].s_id;
        if (tp == 1) {
            x = surf1[id].ts.zent.absi;
            y = surf1[id].ts.zent.ordo;
            z = surf1[id].ts.zent.cote;
            r = surf1[id].ts.rad;
            n1 = msh.entity[i].frvrt;
            dis = wodt_dist_gilq(msh.knot[n1], surf1[id].ts.zent);
            diff = fabs(dis - r);
            err = err + diff;
            nb++;
        }
        if (tp == 2) {
            x = -1.0;
            y = -1.0;
            z = -1.0;
            r = -1.0;
        }
        if (tp == 3) {
            w = surf3[id].pc.sitn;
            x = -1.0;
            y = -1.0;
            z = -1.0;
            r = -1.0;
        }
        fprintf(fp, "%d  %d  %d  %f  %f  %f  %f  %.16f  %.16f  %.16f\n", msh.entity[i].frvrt, msh.entity[i].scvrt, msh.entity[i].thvrt, col[i].red, col[i].green, col[i].blue, x, y, z, r);
    }
    fclose(fp);
    err = err / (double) nb;
    fprintf(tmpout, "Average error=%e\n", err);
}


void qazf_popu_fawj(manif_tl msh, rgb_lk * col, trmsrf * surf1, trmsrf * surf2, trmsrf * surf3, supp_surf * ssr, manif_tl * clmns, rgb_lk * col_f, sphere * S)
{
    int nel, i, tp, id, w, n1, nb;
    double x, y, z, r, dis, diff, err;
    nel = msh.e_grs;
    err = 0.0;
    nb = 0;
    for (i = 0; i < nel; i++) {
        tp = ssr[i].s_type;
        id = ssr[i].s_id;
        if (tp == 1) {
            x = surf1[id].ts.zent.absi;
            y = surf1[id].ts.zent.ordo;
            z = surf1[id].ts.zent.cote;
            r = surf1[id].ts.rad;
            n1 = msh.entity[i].frvrt;
            dis = wodt_dist_gilq(msh.knot[n1], surf1[id].ts.zent);
            diff = fabs(dis - r);
            err = err + diff;
            nb++;
        }
        if (tp == 2) {
            x = -1.0;
            y = -1.0;
            z = -1.0;
            r = -1.0;
        }
        if (tp == 3) {
            w = surf3[id].pc.sitn;
            x = -1.0;
            y = -1.0;
            z = -1.0;
            r = -1.0;
        }
        clmns->entity[i].frvrt = msh.entity[i].frvrt;
        clmns->entity[i].scvrt = msh.entity[i].scvrt;
        clmns->entity[i].thvrt = msh.entity[i].thvrt;
        col_f[i].red = col[i].red;
        col_f[i].green = col[i].green;
        col_f[i].blue = col[i].blue;
        S[i].zent.absi = x;
        S[i].zent.ordo = y;
        S[i].zent.cote = z;
        S[i].rad = r;
    }
    clmns->e_grs = nel;
    err = err / (double) nb;
    fprintf(tmpout, "Average error=%e\n", err);
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

