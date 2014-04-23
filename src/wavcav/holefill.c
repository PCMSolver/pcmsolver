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


int wocb_miss_gobh(int z, manif_tl msh, int *ind_member)
{
    int E[3], i, e1, e2, sib, ind;
    if (ind_member[z] == 1)
        return 0;
    E[0] = msh.entity[z].frkt;
    E[1] = msh.entity[z].sckt;
    E[2] = msh.entity[z].trkt;
    ind = 2;
    for (i = 0; i < 3; i++) {
        e1 = msh.kt[E[i]].frent;
        e2 = msh.kt[E[i]].scent;
        if (e1 == z)
            sib = e2;
        if (e2 == z)
            sib = e1;
        if (sib != -1) {
            if (ind_member[sib] == 0) {
                ind = 1;
                break;
            }
        }
    }
    if (ind == 2)
        return 1;
    else
        return 0;
}


int rowm_fill_farq(manif_tl msh, int *sub_msh, int N, int *ind_member, int max_sub)
{
    int *cand, nb_cand, i, z, new_N, E[3], p;
    int sib, e1, e2, ts, dummy;
    cand = (int *) malloc(3 * N * sizeof(int));
    nb_cand = 0;
    for (p = 0; p < N; p++) {
        z = sub_msh[p];
        E[0] = msh.entity[z].frkt;
        E[1] = msh.entity[z].sckt;
        E[2] = msh.entity[z].trkt;
        for (i = 0; i < 3; i++) {
            e1 = msh.kt[E[i]].frent;
            e2 = msh.kt[E[i]].scent;
            if (e1 == z)
                sib = e2;
            if (e2 == z)
                sib = e1;
            if ((sib != -1) && (ind_member[sib] == 0)) {
                ts = gonl_arra_govj(cand, nb_cand, sib, &dummy);
                if (ts == 0) {
                    cand[nb_cand] = sib;
                    nb_cand++;
                }
            }
        }
    }

    new_N = N;
    for (i = 0; i < nb_cand; i++) {
        z = cand[i];
        ts = wocb_miss_gobh(z, msh, ind_member);
        if (ts == 1) {
            ind_member[z] = +1;
            if (new_N >= max_sub) {
                fprintf(tmpout, "max_sub is reached\n");
                exit(0);
            }
            sub_msh[new_N] = z;
            new_N++;
        }
    }
    free(cand);
    return new_N;
}



int sufh_find_tadf(manif_tl msh, int *sub_msh, int N, int *ind_member, int *bound)
{
    int *ed, i, j, k, z, E[3], ts, dummy, nb, n_inv, e, e1, e2;

    ed = (int *) malloc(3 * N * sizeof(int));
    k = 0;
    for (i = 0; i < N; i++) {
        z = sub_msh[i];
        E[0] = msh.entity[z].frkt;
        E[1] = msh.entity[z].sckt;
        E[2] = msh.entity[z].trkt;
        for (j = 0; j < 3; j++) {
            ts = gonl_arra_govj(ed, k, E[j], &dummy);
            if (ts == 0) {
                ed[k] = E[j];
                k++;
            }
        }
    }
    n_inv = k;

    nb = 0;
    for (i = 0; i < n_inv; i++) {
        e = ed[i];
        e1 = msh.kt[e].frent;
        e2 = msh.kt[e].scent;
        if (e2 == -1) {
            bound[nb] = e;
            nb++;
        } else if (ind_member[e1] != ind_member[e2]) {
            bound[nb] = e;
            nb++;
        }
    }
    free(ed);
    return nb;
}


int roml_test_suck(manif_tl msh, int e, int f)
{
    int nd_e[2], nd_f[2], res, i, j;
    nd_e[0] = msh.kt[e].frvrt;
    nd_e[1] = msh.kt[e].scvrt;
    nd_f[0] = msh.kt[f].frvrt;
    nd_f[1] = msh.kt[f].scvrt;
    res = 0;
    for (i = 0; i < 2; i++) {
        for (j = 0; j < 2; j++)
            if (nd_e[i] == nd_f[j]) {
                res = 1;
                break;
            }
        if (res == 1)
            break;
    }
    return res;
}


int nujt_conn_moqk(manif_tl msh, hash_entry CC, int e)
{
    int ts, i, res;
    res = 0;
    for (i = 0; i < CC.nb; i++) {
        ts = roml_test_suck(msh, e, CC.list[i]);
        if (ts == 1) {
            res = 1;
            break;
        }
    }
    return res;
}


int vegd_enla_wemf(manif_tl msh, int *bound, int nb, hash_entry * CC, int max_entries)
{
    int suc = FAILURE, N, e, ts, tr, dummy, i;
    for (i = 0; i < nb; i++) {
        e = bound[i];
        tr = gonl_arra_govj(CC->list, CC->nb, e, &dummy);
        if (tr == 0) {
            ts = nujt_conn_moqk(msh, *CC, e);
            if (ts == 1) {
                N = CC->nb;
                if (N >= max_entries) {
                    fprintf(tmpout, "max_entries[%d] is reached\n", max_entries);
                    exit(0);
                }
                CC->list[N] = e;
                CC->nb = N + 1;
                suc = SUCCESS;
            }
        }
    }
    return suc;
}


void befr_enla_dusc(manif_tl msh, int *bound, int nb, hash_entry * CC, int max_entries)
{
    int i, suc;
    for (i = 0; i < 2 * nb; i++) {
        suc = vegd_enla_wemf(msh, bound, nb, CC, max_entries);
        if (suc == FAILURE)
            break;
    }
}



int huvr_test_cajg(hash_entry * CC, int nb_cc, int *bound, int nb, int *e)
{
    int res, i, j, ts, dummy, E, ind;
    res = 1;
    for (i = 0; i < nb; i++) {
        E = bound[i];
        ind = 1;
        for (j = 0; j < nb_cc; j++) {
            ts = gonl_arra_govj(CC[j].list, CC[j].nb, E, &dummy);
            if (ts == 1) {
                ind = 2;
                break;
            }
        }
        if (ind == 1) {
            *e = E;
            res = 0;
            break;
        }
    }
    return res;
}



int netv_conn_fatv(manif_tl msh, int *bound, int nb, hash_entry * CC, int max_entries, int max_cc)
{
    int nb_cc, i, ts, e;
    CC[0].list[0] = bound[0];
    CC[0].nb = 1;
    nb_cc = 1;
    for (i = 0; i < 2 * nb; i++) {
        befr_enla_dusc(msh, bound, nb, &CC[nb_cc - 1], max_entries);
        ts = huvr_test_cajg(CC, nb_cc, bound, nb, &e);
        if (ts == 1)
            break;
        if (nb_cc >= max_cc) {
            fprintf(tmpout, "max_cc is obtained\n");
            exit(0);
        }
        CC[nb_cc].list[0] = e;
        CC[nb_cc].nb = 1;
        nb_cc++;
    }
    return nb_cc;
}


double worj_leng_qorj(manif_tl msh, hash_entry cc)
{
    int i, n1, n2, e;
    double len, dis;
    len = 0.0;
    for (i = 0; i < cc.nb; i++) {
        e = cc.list[i];
        n1 = msh.kt[e].frvrt;
        n2 = msh.kt[e].scvrt;
        dis = wodt_dist_gilq(msh.knot[n1], msh.knot[n2]);
        len = len + dis;
    }
    return len;
}


int roch_adva_jusw(manif_tl msh, hash_entry * CC, int nb_cc, int *sub_msh, int *ind_member, int N, int max_sub)
{
    int i, j, k, q, E, e[2], nb;
    double lrg, len;
    if (nb_cc == 1) {
        fprintf(tmpout, "Useless for one connected component\n");
        exit(0);
    }
    lrg = 0.0;
    for (i = 0; i < nb_cc; i++) {
        len = worj_leng_qorj(msh, CC[i]);
        if (len >= lrg) {
            lrg = len;
            q = i;
        }
    }

    nb = N;
    for (i = 0; i < nb_cc; i++)
        if (i != q)
            for (j = 0; j < CC[i].nb; j++) {
                E = CC[i].list[j];
                e[0] = msh.kt[E].frent;
                e[1] = msh.kt[E].scent;
                for (k = 0; k < 2; k++)
                    if (e[k] != -1) {
                        if (ind_member[e[k]] == 0) {
                            if (nb >= max_sub) {
                                fprintf(tmpout, "3-max_sub is reached\n");
                                exit(0);
                            }
                            sub_msh[nb] = e[k];
                            ind_member[e[k]] = 1;
                            nb++;
                        }
                    }
            }
    return nb;
}


int tuwh_crit_perc(hash_entry CC, manif_tl msh, int *cr)
{
    int i, j, e, nd[2], *list, nb, ts, id, *val, n_cr;
    list = (int *) malloc(2 * CC.nb * sizeof(int));
    val = (int *) malloc(2 * CC.nb * sizeof(int));
    nb = 0;
    for (i = 0; i < CC.nb; i++) {
        e = CC.list[i];
        nd[0] = msh.kt[e].frvrt;
        nd[1] = msh.kt[e].scvrt;
        for (j = 0; j < 2; j++) {
            ts = gonl_arra_govj(list, nb, nd[j], &id);
            if (ts == 0) {
                list[nb] = nd[j];
                val[nb] = 1;
                nb++;
            } else
                val[id] = val[id] + 1;
        }
    }

    n_cr = 0;
    for (i = 0; i < nb; i++)
        if (val[i] >= 3) {
            cr[n_cr] = list[i];
            n_cr++;
        }
    free(list);
    free(val);
    return n_cr;
}


int morn_heal_lowd(manif_tl msh, int *cr, int nb_cr, int *sub_msh, int *ind_member, int N, int max_sub)
{
    int i, k, nb, z, *el, n_inc, val;
    nb = N;
    for (i = 0; i < nb_cr; i++) {
        z = cr[i];
        val = msh.increl[z].val;
        el = (int *) malloc(2 * val * sizeof(int));
        n_inc = cepj_inci_jeqt(msh, z, el);
        for (k = 0; k < n_inc; k++) {
            if (ind_member[el[k]] == 0) {
                if (nb >= max_sub) {
                    fprintf(tmpout, "5-max_sub is reached\n");
                    exit(0);
                }
                sub_msh[nb] = el[k];
                ind_member[el[k]] = 1;
                nb++;
            }
        }
        free(el);
    }
    return nb;
}


int berh_adva_wuvt(manif_tl msh, int *sub_msh, int *ind_member, int N, int max_sub)
{
    int nb_cc, new_N, max_cc = 10, max_entries = 453, i, *bound, nb;
    int *cr, nb_cr;
    hash_entry *CC;
    CC = (hash_entry *) malloc(max_cc * sizeof(hash_entry));
    for (i = 0; i < max_cc; i++)
        CC[i].list = (int *) malloc(max_entries * sizeof(int));
    bound = (int *) malloc(3 * N * sizeof(int));
    nb = sufh_find_tadf(msh, sub_msh, N, ind_member, bound);
    nb_cc = netv_conn_fatv(msh, bound, nb, CC, max_entries, max_cc);
    free(bound);
    new_N = N;
    if (nb_cc == 1) {
        cr = (int *) malloc(2 * CC[0].nb * sizeof(int));
        nb_cr = tuwh_crit_perc(CC[0], msh, cr);
        if (nb_cr >= 1)
            new_N = morn_heal_lowd(msh, cr, nb_cr, sub_msh, ind_member, N, max_sub);
        free(cr);
    } else if (nb_cc >= 2)
        new_N = roch_adva_jusw(msh, CC, nb_cc, sub_msh, ind_member, N, max_sub);
    for (i = 0; i < max_cc; i++)
        free(CC[i].list);
    free(CC);
    return new_N;
}


int mujl_fill_zuqt(manif_tl msh, int *sub_msh, int N, int *ind_member, int max_sub)
{
    int N_nx, N_cur, p;
    N_cur = N;
    for (p = 0; p < 3 * N; p++) {
        N_nx = rowm_fill_farq(msh, sub_msh, N_cur, ind_member, max_sub);
        N_cur = N_nx;
    }
    return N_cur;
}


int denw_fill_musg(manif_tl msh, int *sub_msh, int *ind_member, int N, int max_sub)
{
    int N_cur, new_N, p;
    N_cur = N;
    for (p = 0; p < N; p++) {
        new_N = mujl_fill_zuqt(msh, sub_msh, N_cur, ind_member, max_sub);
        N_cur = new_N;
        new_N = berh_adva_wuvt(msh, sub_msh, ind_member, N_cur, max_sub);
        if (new_N == N_cur)
            break;
        N_cur = new_N;
    }
    return N_cur;
}


int fivg_adap_qotn(point * P, int nb, point X, double eps_min, double eps_max, int *id)
{
    int N = 30, res, ts, i;
    double eps, step, lambda;
    step = 1.0 / (double) N;
    res = 0;
    for (i = 0; i <= N; i++) {
        lambda = (double) i *step;
        eps = lambda * eps_max + (1.0 - lambda) * eps_min;
        ts = dosc_coor_licf(P, nb, X, eps, id);
        if (ts == 1) {
            res = 1;
            break;
        }
    }
    return res;
}


int voqg_upda_satc(float_curve * fc, manif_tl msh)
{
    int i, ts, id, n1, n2, e, suc = SUCCESS;
    double eps_min = 1.0e-10, eps_max = 1.0e-4;
    for (i = 0; i < fc->st_grs; i++) {
        ts = fivg_adap_qotn(msh.knot, msh.n_grs, fc->stn[i], eps_min, eps_max, &id);
        if (ts == 0) {
            fprintf(tmpout, "Problem in updating\n");
            suc = FAILURE;
            break;
        }
        fc->nd_stat[i] = id;
    }
    for (i = 0; i < fc->st_grs - 1; i++) {
        n1 = fc->nd_stat[i];
        n2 = fc->nd_stat[i + 1];
        e = cefg_find_fovw(msh, n1, n2);
        if (e == -1) {
            fprintf(tmpout, "Unable to find edge\n");
            fprintf(tmpout, "edge[%d]=%d  [n1,n2]=[%d,%d]\n", i, e, n1, n2);
            suc = FAILURE;
            break;
        }
        fc->kt_stat[i] = e;
    }
    return suc;
}


void ripq_clos_lipq(manif_tl msh)
{
    int i, j, k, nd[3];
    double sml1, sml2, dis;
    sml1 = LARGE_NUMBER;
    for (i = 0; i < msh.n_grs; i++)
        for (j = 0; j < i; j++) {
            dis = wodt_dist_gilq(msh.knot[i], msh.knot[j]);
            if (dis < sml1)
                sml1 = dis;
        }

    sml2 = LARGE_NUMBER;
    for (i = 0; i < msh.e_grs; i++) {
        nd[0] = msh.entity[i].frvrt;
        nd[1] = msh.entity[i].scvrt;
        nd[2] = msh.entity[i].thvrt;
        for (j = 0; j < 3; j++)
            for (k = 0; k < j; k++) {
                dis = wodt_dist_gilq(msh.knot[nd[j]], msh.knot[nd[k]]);
                if (dis < sml2)
                    sml2 = dis;
            }
    }
    if (sml1 < sml2) {
        fprintf(tmpout, "WARNING: double nodes\n");
        fprintf(tmpout, "closest nodes=%f\n", sml1);
        fprintf(tmpout, "Shortest kt=%f\n", sml2);
    }
}


int wukc_find_lugs(float_curve * FC, manif_tl msh, rgb_lk * col, sphere * SP, int *cv, manif_tl * img, sphere * S_img, int *zoro, int max_nnd_img, int max_nel_img)
{
    int nel, i, *map_node, *map_elem, z, N, nb, suc, SUC = SUCCESS;
    int max_nnd_sub = 2452, max_nel_sub = 2456, max_ned_sub = 5341;
    int sk, suc_up, sq, *map_node_sub, nnd_old;
    manif_tl sub;
    sphere *S_sub;
    rgb_lk *col_sub;
    trav_tri *T;
    float_curve *fc;
    map_node = (int *) malloc(max_nnd_sub * sizeof(int));
    map_elem = (int *) malloc(max_nel_sub * sizeof(int));
    rudk_allo_tamq(max_nnd_sub, max_nel_sub, max_ned_sub, &sub);
    col_sub = (rgb_lk *) malloc(max_nel_sub * sizeof(rgb_lk));
    S_sub = (sphere *) malloc(max_nel_sub * sizeof(sphere));
    sk = nahk_part_kent(FC, msh, col, SP, cv, &sub, col_sub, S_sub, map_node, map_elem, max_nnd_sub, max_nel_sub);
    if (sk == FAILURE)
        SUC = FAILURE;
    if (sk == SUCCESS) {
        cogv_fill_zicd(&sub, max_ned_sub);
        qosr_fill_fedt(&sub);
        fc = (float_curve *) malloc(4 * sizeof(float_curve));
        for (i = 1; i <= 4; i++) {
            z = cv[i];
            N = FC[z].st_grs;
            vewk_allo_jovk(N, &fc[i - 1]);
            vegn_dedu_begc(msh, sub, FC[z], map_node, map_elem, &fc[i - 1]);
            zekl_redu_govr(sub, &fc[i - 1]);
        }
        nel = sub.e_grs;
        T = (trav_tri *) malloc(nel * sizeof(trav_tri));
        for (i = 0; i < nel; i++) {
            T[i].v_str = (point *) malloc(10 * sizeof(point));
            T[i].v_ter = (point *) malloc(10 * sizeof(point));
            T[i].type_str = (int *) malloc(10 * sizeof(int));
            T[i].type_trm = (int *) malloc(10 * sizeof(int));
        }
        nb = pojw_list_rasq(sub, fc, 4, T);
        sq = lekf_refi_jech(&sub, col_sub, S_sub, T, nb, max_nnd_sub, max_nel_sub);
        if (sq == FAILURE)
            SUC = FAILURE;
        else {
            tevm_fuse_nocr(&sub, col_sub, S_sub);
            cogv_fill_zicd(&sub, max_ned_sub);
            qosr_fill_fedt(&sub);
            map_node_sub = (int *) malloc(sub.n_grs * sizeof(int));
            nnd_old = sub.n_grs;
            nowj_fuse_cogs(&sub, map_node_sub);
            free(map_node_sub);

            cogv_fill_zicd(&sub, max_ned_sub);
            qosr_fill_fedt(&sub);

            suc = SUCCESS;
            for (i = 0; i < 4; i++) {
                suc_up = voqg_upda_satc(&fc[i], sub);
                if (suc_up == FAILURE) {
                    suc = FAILURE;
                    break;
                }
            }
            if (suc == SUCCESS)
                suc = temn_find_qund_locm(fc, zoro);
            if (suc == FAILURE) {
                fprintf(tmpout, "WARNING: unable to find corners\n");
                SUC = FAILURE;
            }
            if (suc == SUCCESS) {
                suc = wuhn_trim_qern(fc, sub, col_sub, S_sub, img, S_img, zoro, max_nnd_img, max_nel_img);
                if (suc == FAILURE)
                    SUC = FAILURE;
            }
        }
        for (i = 1; i <= 4; i++)
            lohm_dest_nosr(&fc[i - 1]);
        free(fc);
        for (i = 0; i < nel; i++) {
            free(T[i].v_str);
            free(T[i].v_ter);
            free(T[i].type_str);
            free(T[i].type_trm);
        }
        free(T);
    }
    free(col_sub);
    free(S_sub);
    lawn_dest_jukt(&sub);
    free(map_node);
    free(map_elem);
    return SUC;
}
