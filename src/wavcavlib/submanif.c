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
#include "cavity.h"
#include "splinemol.h"
#include "geodesic.h"
#include "pln_sph.h"
#include "sas.h"
#include "smooth.h"



int luct_boun_qifd(float_curve * FC, int *cv, int *trav, int max_trav)
{
    int nb, i, j, z, nb_st;
    nb = 0;
    for (i = 1; i <= 4; i++) {
        z = cv[i];
        nb_st = FC[z].st_grs;
        for (j = 0; j < nb_st - 1; j++) {
            if (nb >= max_trav) {
                fprintf(tmpout, "max_trav[%d] is reached\n", max_trav);
                exit(0);
            }
            trav[nb] = FC[z].trv[j];
            nb++;
        }
    }
    return nb;
}


int sotg_expa_hitc(manif_tl msh, int *sub_msh, int *ind_member, int N, int *trav, int nb_trav, int **exc_ed, int max_sub, int *suc)
{
    int nb, i, j, k, ed[3], s, f, e[2], dummy, ts;
    nb = N;
    *suc = FAILURE;
    for (i = 0; i < N; i++) {
        s = sub_msh[i];
        ed[0] = msh.entity[s].frkt;
        ed[1] = msh.entity[s].sckt;
        ed[2] = msh.entity[s].trkt;
        for (j = 0; j < 3; j++)
            if (exc_ed[i][j] == 0) {
                f = ed[j];
                e[0] = msh.kt[f].frent;
                e[1] = msh.kt[f].scent;
                for (k = 0; k < 2; k++)
                    if ((e[k] != -1) && (e[k] != s)) {
                        ts = gonl_arra_govj(trav, nb_trav, e[k], &dummy);
                        if (ts == 0) {
                            if (ind_member[e[k]] == 0) {
                                if (nb >= max_sub) {
                                    fprintf(tmpout, "1-max_sub=%d is reached\n", max_sub);
                                    exit(0);
                                }
                                exc_ed[i][j] = 1;
                                sub_msh[nb] = e[k];
                                ind_member[e[k]] = 1;
                                nb++;
                                *suc = SUCCESS;
                            }
                        }
                    }
            }
    }
    return nb;
}


int remj_enla_lohp(manif_tl msh, int *sub_msh, int *ind_member, int N, int max_sub, int *suc)
{
    int nb, i, j, k, ed[3], s, f, e[2];
    nb = N;
    *suc = FAILURE;
    for (i = 0; i < N; i++) {
        s = sub_msh[i];
        ed[0] = msh.entity[s].frkt;
        ed[1] = msh.entity[s].sckt;
        ed[2] = msh.entity[s].trkt;
        for (j = 0; j < 3; j++) {
            f = ed[j];
            e[0] = msh.kt[f].frent;
            e[1] = msh.kt[f].scent;
            for (k = 0; k < 2; k++)
                if ((e[k] != -1) && (e[k] != s)) {
                    if (ind_member[e[k]] == 0) {
                        if (nb >= max_sub) {
                            fprintf(tmpout, "3-max_sub is reached\n");
                            exit(0);
                        }
                        sub_msh[nb] = e[k];
                        ind_member[e[k]] = 1;
                        nb++;
                        *suc = SUCCESS;
                    }
                }
        }
    }
    return nb;
}



int rukm_expa_cern(int seed, manif_tl msh, int *trav, int nb_trav, int *sub_msh, int *ind_member, int max_rec, int **exc_ed, int max_sub, int *N)
{
    int nb, j, nb_new, nel, suc, SUC = SUCCESS;
    sub_msh[0] = seed;
    ind_member[seed] = 1;
    nb = 1;
    nel = msh.e_grs;
    for (j = 0; j < nel; j++) {
        nb_new = sotg_expa_hitc(msh, sub_msh, ind_member, nb, trav, nb_trav, exc_ed, max_sub, &suc);
        nb = nb_new;
        if (suc == FAILURE)
            break;
        if (nb >= max_rec) {
            SUC = FAILURE;
            break;
        }
    }
    *N = nb;
    return SUC;
}



void limr_simp_wufg(manif_tl msh, int *trav, int nb_trav, int *sub_msh, int *ind_member, int max_sub, int *NB)
{
    int N, i, suc_dummy, new_N;
    N = *NB;
    for (i = 0; i < nb_trav; i++) {
        if (N >= max_sub) {
            fprintf(tmpout, "2-max_sub is reached: N=%d\n", N);
            exit(0);
        }
        sub_msh[N] = trav[i];
        ind_member[trav[i]] = 1;
        N++;
    }
    new_N = remj_enla_lohp(msh, sub_msh, ind_member, N, max_sub, &suc_dummy);
    N = new_N;
    *NB = N;
}



int tedk_mail_wetr(manif_tl msh, int *trav, int nb_trav, int *sub_msh, int *ind_member, int max_sub, int *NB)
{
    int max_rec, nel, *seed, i, j, k, nb_sd, suc, *exc_seed;
    int ed[3], e[2], ts, tr, dummy, s, f, ind, N, **exc_ed, p, q;
    nel = msh.e_grs;
    max_rec = 1000;

    seed = (int *) malloc(3 * nb_trav * sizeof(int));
    nb_sd = 0;
    for (i = 0; i < nb_trav; i++) {
        s = trav[i];
        ed[0] = msh.entity[s].frkt;
        ed[1] = msh.entity[s].sckt;
        ed[2] = msh.entity[s].trkt;
        for (j = 0; j < 3; j++) {
            f = ed[j];
            e[0] = msh.kt[f].frent;
            e[1] = msh.kt[f].scent;
            for (k = 0; k < 2; k++)
                if (e[k] != -1) {
                    ts = gonl_arra_govj(trav, nb_trav, e[k], &dummy);
                    if (ts == 0) {
                        tr = gonl_arra_govj(seed, nb_sd, e[k], &dummy);
                        if (tr == 0) {
                            seed[nb_sd] = e[k];
                            nb_sd++;
                        }
                    }
                }
        }
    }

    exc_seed = (int *) malloc(nb_sd * sizeof(int));
    for (i = 0; i < nb_sd; i++)
        exc_seed[i] = 0;
    exc_ed = (int **) malloc(max_sub * sizeof(int *));
    for (i = 0; i < max_sub; i++)
        exc_ed[i] = (int *) malloc(3 * sizeof(int));
    ind = 1;
    for (i = 0; i < nb_sd; i++)
        if (exc_seed[i] == 0) {
            for (p = 0; p < max_sub; p++)
                for (q = 0; q < 3; q++)
                    exc_ed[p][q] = 0;
            for (p = 0; p < msh.e_grs; p++)
                ind_member[p] = 0;
            suc = rukm_expa_cern(seed[i], msh, trav, nb_trav, sub_msh, ind_member, max_rec, exc_ed, max_sub, &N);
            if (suc == SUCCESS) {
                ind = 2;
                break;
            }
            for (j = 0; j < nb_sd; j++)
                if (exc_seed[j] == 0) {
                    if (ind_member[seed[j]] == 1)
                        exc_seed[j] = 1;
                }
        }
    free(seed);
    free(exc_seed);
    for (i = 0; i < max_sub; i++)
        free(exc_ed[i]);
    free(exc_ed);
    if (ind == 1)
        return FAILURE;
    *NB = N;
    return SUCCESS;
}


int husm_elem_gumn(float_curve * FC, int *cv, manif_tl msh, int *sub_msh, int *ind_member, int max_sub, int *NB)
{
    int N, *trav, nb_trav, max_trav = 405, i, ts;
    int dummy, suc, suc_dummy, new_N;
    trav = (int *) malloc(max_trav * sizeof(int));
    nb_trav = luct_boun_qifd(FC, cv, trav, max_trav);
    suc = tedk_mail_wetr(msh, trav, nb_trav, sub_msh, ind_member, max_sub, &N);
    if (suc == SUCCESS) {
        for (i = 0; i < nb_trav; i++) {
            ts = gonl_arra_govj(sub_msh, N, trav[i], &dummy);
            if (ts == 0) {
                if (N >= max_sub) {
                    fprintf(tmpout, "2-max_sub is reached: N=%d\n", N);
                    exit(0);
                }
                sub_msh[N] = trav[i];
                ind_member[trav[i]] = 1;
                N++;
            }
        }
    }
    if (suc == FAILURE) {
        N = 0;
        for (i = 0; i < msh.e_grs; i++)
            ind_member[i] = 0;
        limr_simp_wufg(msh, trav, nb_trav, sub_msh, ind_member, max_sub, &N);
        suc = SUCCESS;
    }

    new_N = remj_enla_lohp(msh, sub_msh, ind_member, N, max_sub, &suc_dummy);
    N = new_N;

    new_N = denw_fill_musg(msh, sub_msh, ind_member, N, max_sub);
    N = new_N;
    free(trav);
    *NB = N;
    return suc;
}


int gepj_comm_jasq(float_curve FC1, float_curve FC2, int *cm)
{
    int nd1[2], nd2[2], n, res, i, j;
    n = FC1.st_grs;
    nd1[0] = FC1.nd_idx[0];
    nd1[1] = FC1.nd_idx[n - 1];

    n = FC2.st_grs;
    nd2[0] = FC2.nd_idx[0];
    nd2[1] = FC2.nd_idx[n - 1];

    res = 0;
    for (i = 0; i < 2; i++) {
        for (j = 0; j < 2; j++)
            if (nd1[i] == nd2[j]) {
                res = 1;
                *cm = nd1[i];
                break;
            }
        if (res == 1)
            break;
    }
    return res;
}



int temn_find_qund_locm(float_curve * fc, int *zoro)
{
    int suc = FAILURE, i, j, nb, ts, cm;
    nb = 0;
    for (i = 0; i < 4; i++)
        for (j = 0; j < i; j++) {
            ts = gepj_comm_jasq(fc[i], fc[j], &cm);
            if (ts == 1) {
                zoro[nb] = cm;
                nb++;
            }
        }
    if (nb == 4)
        suc = SUCCESS;
    return suc;
}


int nahk_part_kent(float_curve * FC, manif_tl msh, rgb_lk * col, sphere * S, int *cv, manif_tl * sub, rgb_lk * col_sub, sphere * S_sub, int *map_node, int *map_elem, int max_node_sub, int max_elem_sub)
{
    int *sub_msh, nb_sub, i, max_sub;
    int z, suc, nel, *ind_member;
    nel = msh.e_grs;

    max_sub = nel / 2;
    sub_msh = (int *) malloc(max_sub * sizeof(int));
    ind_member = (int *) malloc(nel * sizeof(int));
    suc = husm_elem_gumn(FC, cv, msh, sub_msh, ind_member, max_sub, &nb_sub);
    if (suc == SUCCESS) {
        tugw_extr_vuzf(msh, sub, sub_msh, nb_sub, map_node, map_elem, max_node_sub, max_elem_sub);
        for (i = 0; i < sub->e_grs; i++) {
            z = map_elem[i];
            col_sub[i].red = col[z].red;
            col_sub[i].green = col[z].green;
            col_sub[i].blue = col[z].blue;
            neqg_find_lodr_bogm(S[z], &S_sub[i]);
        }
    }

    free(sub_msh);
    free(ind_member);

    return suc;
}


void vafh_corr_monz(manif_tl sub, int i, float_curve fc_in, float_curve * fc_out, int *map_node)
{
    int nd, ts, id;
    nd = fc_in.nd_idx[i];
    ts = gonl_arra_govj(map_node, sub.n_grs, nd, &id);
    if (ts == 0) {
        fprintf(tmpout, "Unable to associate parent node\n");
        exit(0);
    }
    fc_out->nd_idx[i] = id;
}


void qenc_corr_tunq(manif_tl msh, manif_tl sub, int i, float_curve fc_in, float_curve * fc_out, int *map_node)
{
    int n1, n2, e, e_sub, ts1, ts2, id1, id2;
    e = fc_in.kt_idx[i];
    n1 = msh.kt[e].frvrt;
    n2 = msh.kt[e].scvrt;
    ts1 = gonl_arra_govj(map_node, sub.n_grs, n1, &id1);
    ts2 = gonl_arra_govj(map_node, sub.n_grs, n2, &id2);
    if ((ts1 == 0) || (ts2 == 0)) {
        fprintf(tmpout, "Unable to associate parent edge\n");
        exit(0);
    }
    e_sub = cefg_find_fovw(sub, id1, id2);
    fc_out->kt_idx[i] = e_sub;
}


void vegn_dedu_begc(manif_tl msh, manif_tl sub, float_curve fc_in, int *map_node, int *map_elem, float_curve * fc_out)
{
    int N, i, cs, z, ts, id;
    N = fc_in.st_grs;
    for (i = 0; i < N; i++)
        getf_find_rogc_todj(fc_in.stn[i], &fc_out->stn[i]);
    for (i = 0; i < N - 1; i++) {
        cs = fc_in.cs[i];
        fc_out->cs[i] = cs;
        if ((cs == 1) || (cs == 2)) {
            vafh_corr_monz(sub, i, fc_in, fc_out, map_node);
            if ((i == N - 2) && (cs == 1))
                vafh_corr_monz(sub, N - 1, fc_in, fc_out, map_node);
            if ((i == N - 2) && (cs == 2))
                qenc_corr_tunq(msh, sub, N - 1, fc_in, fc_out, map_node);
        }
        if ((cs == 3) || (cs == 4)) {
            qenc_corr_tunq(msh, sub, i, fc_in, fc_out, map_node);
            if ((i == N - 2) && (cs == 3))
                qenc_corr_tunq(msh, sub, N - 1, fc_in, fc_out, map_node);
            if ((i == N - 2) && (cs == 4))
                vafh_corr_monz(sub, N - 1, fc_in, fc_out, map_node);
        }
    }

    for (i = 0; i < N - 1; i++) {
        z = fc_in.trv[i];
        ts = gonl_arra_govj(map_elem, sub.e_grs, z, &id);
        fc_out->trv[i] = id;
    }
    fc_out->st_grs = N;
}


int topr_exte_jemt(manif_tl msh, int *sub_msh, int N, int *fence, int nb_fence, int *suc)
{
    int nb, i, j, k, ed[4], s, f, ts, tr, e[2], dummy;
    nb = N;
    *suc = FAILURE;
    for (i = 0; i < N; i++) {
        s = sub_msh[i];
        ed[1] = msh.entity[s].frkt;
        ed[2] = msh.entity[s].sckt;
        ed[3] = msh.entity[s].trkt;
        for (j = 1; j <= 3; j++) {
            f = ed[j];
            ts = gonl_arra_govj(fence, nb_fence, f, &dummy);
            if (ts == 0) {
                e[0] = msh.kt[f].frent;
                e[1] = msh.kt[f].scent;
                for (k = 0; k < 2; k++)
                    if ((e[k] != -1) && (e[k] != s)) {
                        tr = gonl_arra_govj(sub_msh, nb, e[k], &dummy);
                        if (tr == 0) {
                            sub_msh[nb] = e[k];
                            nb++;
                            *suc = SUCCESS;
                        }
                    }
            }
        }
    }
    return nb;
}


int csm = 0;


int mudp_corr_rost(int *sub_msh, int nb, manif_tl msh, int *fence, int nb_fence)
{
    int *map_node, *map_elem, nel, n_bound, i, res;
    int n1, n2, *im_1, *im_2, n_a, n_b, z, ind, j, *exc;
    manif_tl temp;
    nel = nb;
    map_elem = (int *) malloc(nel * sizeof(int));
    map_node = (int *) malloc(3 * nel * sizeof(int));
    juqr_allo_nohd(3 * nel, nel, 3 * nel, &temp);
    tugw_extr_vuzf(msh, &temp, sub_msh, nb, map_node, map_elem, 3 * nel, nel);
    cogv_fill_zicd(&temp, 3 * nel);
    im_1 = (int *) malloc(temp.k_grs * sizeof(int));
    im_2 = (int *) malloc(temp.k_grs * sizeof(int));
    n_bound = 0;
    for (i = 0; i < temp.k_grs; i++)
        if (temp.kt[i].scent == -1) {
            n1 = temp.kt[i].frvrt;
            n2 = temp.kt[i].scvrt;
            im_1[n_bound] = map_node[n1];
            im_2[n_bound] = map_node[n2];
            n_bound++;
        }
    free(map_elem);
    free(map_node);
    zesg_dest_cokd(&temp);
    if (n_bound != nb_fence) {
        free(im_1);
        free(im_2);
        return 0;
    } else {
        exc = (int *) malloc(n_bound * sizeof(int));
        for (i = 0; i < n_bound; i++)
            exc[i] = 0;
        res = 1;
        for (i = 0; i < n_bound; i++) {
            n_a = im_1[i];
            n_b = im_2[i];
            ind = 0;
            for (j = 0; j < nb_fence; j++)
                if (exc[j] == 0) {
                    z = fence[j];
                    n1 = msh.kt[z].frvrt;
                    n2 = msh.kt[z].scvrt;
                    if (((n_a == n1) && (n_b == n2)) || ((n_a == n2) && (n_b == n1))) {
                        ind = 1;
                        exc[j] = 1;
                        break;
                    }
                }
            if (ind == 0) {
                res = 0;
                break;
            }
        }
        free(exc);
    }
    free(im_1);
    free(im_2);
    return res;
}


int ticq_trim_pehf(float_curve * fc, manif_tl msh, rgb_lk * col_msh, int *sub_msh, int *NB)
{
    int *fence, nb_fence, i, j, k, N, M, nb, suc, tr, dummy;
    int nb_new, src[2], ind, ts, SUC = SUCCESS;
    M = 0;
    for (i = 0; i < 4; i++)
        M = M + fc[i].st_grs;
    fence = (int *) malloc(M * sizeof(int));

    k = 0;
    for (i = 0; i < 4; i++) {
        N = fc[i].st_grs;
        for (j = 0; j < N - 1; j++) {
            tr = gonl_arra_govj(fence, k, fc[i].kt_stat[j], &dummy);
            if (tr == 0) {
                fence[k] = fc[i].kt_stat[j];
                k++;
            }
        }
    }
    nb_fence = k;

    src[0] = msh.kt[fence[0]].frent;
    src[1] = msh.kt[fence[0]].scent;
    ind = 1;
    for (i = 0; i < 2; i++) {
        if (src[i] != -1) {
            sub_msh[0] = src[i];
            nb = 1;
            while (1) {
                nb_new = topr_exte_jemt(msh, sub_msh, nb, fence, nb_fence, &suc);
                nb = nb_new;
                if (suc == FAILURE)
                    break;
            }
            ts = mudp_corr_rost(sub_msh, nb, msh, fence, nb_fence);
            if (ts == 1) {
                ind = 2;
                break;
            }
        }
    }
    free(fence);
    if (ind == 1) {
        fprintf(tmpout, "Unable to find patch\n");
        SUC = FAILURE;

    }
    *NB = nb;
    return SUC;
}



int wuhn_trim_qern(float_curve * fc, manif_tl msh, rgb_lk * col_msh, sphere * S_msh, manif_tl * img, sphere * S_img, int *zoro, int max_nnd_img, int max_nel_img)
{
    int *sub_msh, nb, nel, *map_node, *map_elem, i, z, ts, id, suc, w;
    nel = msh.e_grs;
    sub_msh = (int *) malloc(nel * sizeof(int));
    suc = ticq_trim_pehf(fc, msh, col_msh, sub_msh, &nb);
    if (suc == SUCCESS) {
        map_elem = (int *) malloc(nel * sizeof(int));
        map_node = (int *) malloc(3 * nel * sizeof(int));
        tugw_extr_vuzf(msh, img, sub_msh, nb, map_node, map_elem, max_nnd_img, max_nel_img);
        for (i = 0; i < 4; i++) {
            z = zoro[i];
            ts = gonl_arra_govj(map_node, img->n_grs, z, &id);
            if (ts == 0) {
                fprintf(tmpout, "Unable to find patch corners\n");

                suc = FAILURE;
                break;
            }
            zoro[i] = id;
        }
        for (i = 0; i < img->e_grs; i++) {
            w = map_elem[i];
            neqg_find_lodr_bogm(S_msh[w], &S_img[i]);
        }
        free(map_node);
        free(map_elem);

    }
    free(sub_msh);
    return suc;
}
