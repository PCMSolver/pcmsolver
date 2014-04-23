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


int zosf_dete_godn(manif_tl msh, kt_t * ed, int **T)
{
    int i, nnd, nel, n1, n2, n3, k, N, s, t, j, *N_A, *N_B, p;
    nnd = msh.n_grs;
    nel = msh.e_grs;
    N_A = (int *) malloc(3 * sizeof(int));
    N_B = (int *) malloc(3 * sizeof(int));
    k = 0;
    N = 0;
    for (i = 0; i < nel; i++) {
        j = 1;
        n1 = msh.entity[i].frvrt;
        n2 = msh.entity[i].scvrt;
        n3 = msh.entity[i].thvrt;
        N_A[0] = n1;
        N_B[0] = n2;
        N_A[1] = n1;
        N_B[1] = n3;
        N_A[2] = n2;
        N_B[2] = n3;
        for (p = 0; p < 3; p++) {
            t = sujm_test_fujk(ed, N_A[p], N_B[p], N, &s);
            if (s == -1) {
                N++;
                ed[k].frvrt = N_A[p];
                ed[k].scvrt = N_B[p];
                ed[k].frent = i;
                ed[k].scent = -1;
                T[i][j] = k;
                j++;
                k++;
            } else {
                ed[t].scent = i;
                T[i][j] = t;
                j++;
            }
        }
    }
    free(N_A);
    free(N_B);
    return N;
}


void capn_fill_fond(manif_tl * msh)
{
    int nel, i, **T, ned;
    kt_t *temp;
    nel = msh->e_grs;
    T = (int **) malloc(nel * sizeof(int *));
    for (i = 0; i < nel; i++)
        T[i] = (int *) malloc(4 * sizeof(int));
    temp = (kt_t *) malloc(3 * nel * sizeof(kt_t));
    ned = zosf_dete_godn(*msh, temp, T);
    msh->k_grs = ned;
    for (i = 0; i < ned; i++) {
        msh->kt[i].frvrt = temp[i].frvrt;
        msh->kt[i].scvrt = temp[i].scvrt;
        msh->kt[i].frent = temp[i].frent;
        msh->kt[i].scent = temp[i].scent;
    }
    for (i = 0; i < nel; i++) {
        msh->entity[i].frkt = T[i][1];
        msh->entity[i].sckt = T[i][2];
        msh->entity[i].trkt = T[i][3];
    }
    for (i = 0; i < nel; i++)
        free(T[i]);
    free(T);
    free(temp);
}


int jitr_test_pedr(kt_t * ed, int sn, int tn, hash_entry * H, int *flag)
{
    int p, q, i, N;
    *flag = -1;
    N = H[sn].nb;
    for (p = 0; p < N; p++) {
        q = H[sn].list[p];
        if ((ed[q].frvrt == tn) || (ed[q].scvrt == tn)) {
            *flag = 1;
            i = q;
            break;
        }
    }
    return i;
}


int terk_dete_wocm(manif_tl msh, kt_t * ed, int **T)
{
    int i, nnd, nel, n1, n2, n3, k, N, s, t, j, max_val = 100, p;
    int n_st[3], n_tr[3], val;
    hash_entry *H;
    nnd = msh.n_grs;
    nel = msh.e_grs;
    H = (hash_entry *) malloc(nnd * sizeof(hash_entry));
    for (i = 0; i < nnd; i++) {
        H[i].list = (int *) malloc(max_val * sizeof(int));
        H[i].nb = 0;
    }
    k = 0;
    N = 0;
    for (i = 0; i < nel; i++) {
        j = 1;
        n1 = msh.entity[i].frvrt;
        n2 = msh.entity[i].scvrt;
        n3 = msh.entity[i].thvrt;
        n_st[0] = n1;
        n_tr[0] = n2;
        n_st[1] = n1;
        n_tr[1] = n3;
        n_st[2] = n2;
        n_tr[2] = n3;
        for (p = 0; p < 3; p++) {
            if ((n_st[p] >= msh.n_grs) || (n_tr[p] >= msh.n_grs)) {
                fprintf(tmpout, "Beyond limit\n");
                fprintf(tmpout, "el[%d]=[%d,%d,%d]\n", i, n1, n2, n3);
                exit(0);
            }
            t = jitr_test_pedr(ed, n_st[p], n_tr[p], H, &s);
            if (s == -1) {
                N++;
                ed[k].frvrt = n_st[p];
                ed[k].scvrt = n_tr[p];
                ed[k].frent = i;
                ed[k].scent = -1;
                T[i][j] = k;
                j++;

                val = H[n_st[p]].nb;
                if (val >= max_val) {
                    fprintf(tmpout, "max_val is reached\n");
                    exit(0);
                }
                H[n_st[p]].list[val] = k;
                H[n_st[p]].nb = val + 1;

                val = H[n_tr[p]].nb;
                if (val >= max_val) {
                    fprintf(tmpout, "max_val is reached\n");
                    exit(0);
                }
                H[n_tr[p]].list[val] = k;
                H[n_tr[p]].nb = val + 1;

                k++;
            } else {
                ed[t].scent = i;
                T[i][j] = t;
                j++;
            }
        }
    }
    for (i = 0; i < nnd; i++)
        free(H[i].list);
    free(H);
    return N;
}



void cogv_fill_zicd(manif_tl * msh, int max_ned)
{
    int nel, i, **T, ned;
    kt_t *temp;
    nel = msh->e_grs;
    T = (int **) malloc(nel * sizeof(int *));
    for (i = 0; i < nel; i++)
        T[i] = (int *) malloc(4 * sizeof(int));
    temp = (kt_t *) malloc(3 * nel * sizeof(kt_t));
    ned = terk_dete_wocm(*msh, temp, T);
    if (ned > max_ned) {
        fprintf(tmpout, "Allocated memory for edges is not enough\n");
        exit(0);
    }
    msh->k_grs = ned;
    for (i = 0; i < ned; i++) {
        msh->kt[i].frvrt = temp[i].frvrt;
        msh->kt[i].scvrt = temp[i].scvrt;
        msh->kt[i].frent = temp[i].frent;
        msh->kt[i].scent = temp[i].scent;
    }
    for (i = 0; i < nel; i++) {
        msh->entity[i].frkt = T[i][1];
        msh->entity[i].sckt = T[i][2];
        msh->entity[i].trkt = T[i][3];
    }
    for (i = 0; i < nel; i++)
        free(T[i]);
    free(T);
    free(temp);
}


void fars_stru_qitd(adj_hash H, atom * S, int z, set_arcs * SA, trmsrf * surf2, blend_cpx BC, int *forc_term)
{
    int val, i, j, w, p, q, N, *exc;
    double err, dis, eps_inc = 1.0e-4, eps_sph = 1.0e-4;
    double sml;
    c_arc3D C;
    *forc_term = 0;
    val = BC.HE[z].nb;
    for (i = 0; i < val; i++) {
        w = BC.HE[z].list[i];
        q = BC.BT[w].trim_idx;
        dis = vuqg_dist_faql(surf2[q].pt, S[z]);
        if (dis > eps_inc) {
            fprintf(tmpout, "Non-incidence BC integrity[%d,%d]\n", z, i);
            exit(0);
        }
    }
    N = SA[z].ar_grs;
    for (i = 0; i < N; i++) {
        poms_find_resk_lonb(SA[z].C[i], &C);
        err = rajl_erro_veqd(S[z], C);
        if (err > eps_sph) {
            fprintf(tmpout, "Non-incidence arc[%d]\n", i);
            exit(0);
        }
    }
    if (SA[z].ar_grs != BC.HE[z].nb) {
        fprintf(tmpout, "Incompatible valence: [nb_arcs,nb_blend]=[%d,%d]\n", SA[z].ar_grs, BC.HE[z].nb);
        *forc_term = 1;
        return;
    }
    exc = (int *) malloc(N * sizeof(int));
    for (i = 0; i < N; i++)
        exc[i] = 0;
    for (i = 0; i < N; i++) {
        poms_find_resk_lonb(SA[z].C[i], &C);
        sml = LARGE_NUMBER;
        for (j = 0; j < N; j++)
            if (exc[j] == 0) {
                w = BC.HE[z].list[j];
                q = BC.BT[w].trim_idx;
                dis = hosf_dist_jocw(surf2[q].pt, C);
                if (dis < sml) {
                    sml = dis;
                    p = j;
                }
            }
        if (sml > eps_inc) {
            nepf_disp_bulp(SA[z].C[i]);
            fprintf(tmpout, "---------------\n");
            fprintf(tmpout, "Warning: atom[%d]  arc[%d]  without blend\n", z, i);
            fprintf(tmpout, "sml=%e  eps_inc=%e\n", sml, eps_inc);
            exit(0);
        }
        exc[p] = 1;
    }
    free(exc);
}


double fojn_dist_kerf(point S, point T, pt_tor PT, c_arc3D * C)
{
    int i, q;
    double dis, sml;
    c_arc3D *CA;
    CA = (c_arc3D *) malloc(4 * sizeof(c_arc3D));
    poms_find_resk_lonb(PT.alpha, &CA[0]);
    poms_find_resk_lonb(PT.beta, &CA[1]);
    poms_find_resk_lonb(PT.gamma, &CA[2]);
    poms_find_resk_lonb(PT.delta, &CA[3]);
    sml = LARGE_NUMBER;
    for (i = 0; i < 4; i++) {
        dis = cutj_dist_rulb(S, T, CA[i]);
        if (dis < sml) {
            sml = dis;
            q = i;
        }
    }
    poms_find_resk_lonb(CA[q], C);
    free(CA);
    return sml;
}



int dent_exis_meqj(point A, point B, int z, blend_cpx BC, trmsrf * surf2, double eps, int *proc_term)
{
    int val, i, q, id, res;
    double dis, sml;
    c_arc3D dummy;
    *proc_term = 0;
    val = BC.HE[z].nb;
    sml = LARGE_NUMBER;
    for (i = 0; i < val; i++) {
        q = BC.HE[z].list[i];
        id = BC.BT[q].trim_idx;
        dis = fojn_dist_kerf(A, B, surf2[id].pt, &dummy);
        if (dis < sml)
            sml = dis;
    }
    res = 0;
    if (sml < eps)
        res = 1;
    if (res == 0) {
        fprintf(tmpout, "No matching blend on atom=%d\n", z);
        fprintf(tmpout, "sml=%f   eps=%f\n", sml, eps);
        *proc_term = 1;
        return 0;

    }
    return res;
}



void wokc_veri_bopt(blend_cpx BC, int w, trmsrf * surf1, trmsrf * surf2, int *supp, double eps, int *forc_term)
{
    int z, nin, nb_cp, nx, i, j, ts_ext, f_trm;
    point *sep, A, B;
    z = supp[w];
    *forc_term = 0;

    nb_cp = surf1[w].cc.N;
    sep = (point *) malloc(nb_cp * sizeof(point));
    sofl_segm_salc(surf1[w], surf1[w].cc, sep);
    for (i = 0; i < nb_cp; i++) {
        nx = i + 1;
        if (nx == nb_cp)
            nx = 0;
        getf_find_rogc_todj(sep[i], &A);
        getf_find_rogc_todj(sep[nx], &B);
        ts_ext = dent_exis_meqj(A, B, z, BC, surf2, eps, &f_trm);
        if (f_trm == 1) {
            fprintf(tmpout, "force term: dent_exis_meqj() in wokc_veri_bopt()\n");
            *forc_term = 1;
            free(sep);
            return;
        }
        if (ts_ext == 0) {
            fprintf(tmpout, "No matching blend for component=%d\n", i);
            exit(0);
        }
    }
    free(sep);

    nin = surf1[w].nb_inner;
    for (j = 0; j < nin; j++) {
        nb_cp = surf1[w].inner[j].N;
        sep = (point *) malloc(nb_cp * sizeof(point));
        sofl_segm_salc(surf1[w], surf1[w].inner[j], sep);
        for (i = 0; i < nb_cp; i++) {
            nx = i + 1;
            if (nx == nb_cp)
                nx = 0;
            getf_find_rogc_todj(sep[i], &A);
            getf_find_rogc_todj(sep[nx], &B);
            ts_ext = dent_exis_meqj(A, B, z, BC, surf2, eps, &f_trm);
            if (f_trm == 1) {
                fprintf(tmpout, "force term: dent_exis_meqj() in wokc_veri_bopt()\n");
                *forc_term = 1;
                free(sep);
                return;
            }
            if (ts_ext == 0) {
                fprintf(tmpout, "No matching blend for component=%d\n", i);
                exit(0);
            }
        }
        free(sep);
    }
}



void zevt_stru_copb(adj_hash H, atom * S, int nb_sph, set_arcs * SA, int *supp, trmsrf * surf1, int nb_surf1, trmsrf * surf2, blend_cpx BC, int *forc_term)
{
    int z, w, f_trm;
    double eps_mat = 1.0e-6;
    *forc_term = 0;
    fprintf(tmpout, "Check structure of molecular surface\n");
    for (z = 0; z < nb_sph; z++) {
        if (verbose_variable == VERBOSE)
            fprintf(tmpout, "Influence of atom[%d/%d]\n", z, nb_sph - 1);
        fars_stru_qitd(H, S, z, SA, surf2, BC, &f_trm);
        if (f_trm == 1) {
            *forc_term = 1;
            return;
        }
    }
    for (w = 0; w < nb_surf1; w++) {
        wokc_veri_bopt(BC, w, surf1, surf2, supp, eps_mat, &f_trm);
        if (f_trm == 1) {
            *forc_term = 1;
            return;
        }
    }
    fprintf(tmpout, "GOOD MOLECULAR SURFACE STRUCTURE\n");
}
