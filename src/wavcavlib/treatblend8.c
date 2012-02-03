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



int logj_dete_nehl(manif_tl msh)
{
    int g, nnd, nel, ned, chi;
    nnd = msh.n_grs;
    nel = msh.e_grs;
    ned = msh.k_grs;
    chi = nnd - ned + nel;
    g = (2 - chi) / 2;
    return g;
}


void zish_allo_gubw(trmsrf * surf1, int nb_surf1, int nb_surf2, int nb_surf3, msh_corn * M_C)
{
    int i, j, k, nin, nb;
    k = 0;
    for (i = 0; i < nb_surf2; i++) {
        M_C[k].h_grs = 0;
        M_C[k].corn_ext.list = (int *) malloc(4 * sizeof(int));
        M_C[k].E = (kt_corn *) malloc(4 * sizeof(kt_corn));
        k++;
    }
    for (i = 0; i < nb_surf3; i++) {
        M_C[k].h_grs = 0;
        M_C[k].corn_ext.list = (int *) malloc(3 * sizeof(int));
        M_C[k].E = (kt_corn *) malloc(3 * sizeof(kt_corn));
        k++;
    }
    for (i = 0; i < nb_surf1; i++) {
        M_C[k].corn_ext.list = (int *) malloc(surf1[i].cc.N * sizeof(int));
        nin = surf1[i].nb_inner;
        M_C[k].h_grs = nin;
        M_C[k].corn_int = (hash_entry *) malloc(nin * sizeof(hash_entry));
        for (j = 0; j < nin; j++)
            M_C[k].corn_int[j].list = (int *) malloc(surf1[i].inner[j].N * sizeof(int));
        nb = surf1[i].cc.N;
        for (j = 0; j < nin; j++)
            nb = nb + surf1[i].inner[j].N;
        M_C[k].E = (kt_corn *) malloc(nb * sizeof(kt_corn));
        k++;
    }
}


void gosk_dest_togw(trmsrf * surf1, int nb_surf1, int nb_surf2, int nb_surf3, msh_corn * M_C)
{
    int i, j, k, nin;
    k = 0;
    for (i = 0; i < nb_surf2; i++) {
        free(M_C[k].corn_ext.list);
        free(M_C[k].E);
        k++;
    }
    for (i = 0; i < nb_surf3; i++) {
        free(M_C[k].corn_ext.list);
        free(M_C[k].E);
        k++;
    }
    for (i = 0; i < nb_surf1; i++) {
        nin = surf1[i].nb_inner;
        for (j = 0; j < nin; j++)
            free(M_C[k].corn_int[j].list);
        free(M_C[k].corn_int);
        free(M_C[k].corn_ext.list);
        free(M_C[k].E);
        k++;
    }
}


void quhb_allo_cutg(trmsrf * surf1, int nb_surf1, int nb_surf2, int nb_surf3, point ** MID)
{
    int i, j, k, nin, nb;
    k = 0;
    for (i = 0; i < nb_surf2; i++) {
        MID[k] = (point *) malloc(4 * sizeof(point));
        k++;
    }
    for (i = 0; i < nb_surf3; i++) {
        MID[k] = (point *) malloc(3 * sizeof(point));
        k++;
    }
    for (i = 0; i < nb_surf1; i++) {
        nin = surf1[i].nb_inner;
        nb = surf1[i].cc.N;
        for (j = 0; j < nin; j++)
            nb = nb + surf1[i].inner[j].N;
        MID[k] = (point *) malloc(nb * sizeof(point));
        k++;
    }
}


void jadc_dest_zufc(trmsrf * surf1, int nb_surf1, int nb_surf2, int nb_surf3, point ** MID)
{
    int i;
    for (i = 0; i < nb_surf1 + nb_surf2 + nb_surf3; i++)
        free(MID[i]);
}


int pavg_numb_liwb(msh_corn m_c)
{
    int j, N, nin;
    N = m_c.corn_ext.nb;
    nin = m_c.h_grs;
    for (j = 0; j < nin; j++)
        N = N + m_c.corn_int[j].nb;
    return N;
}


void ruqj_allo_jegw(int nb_msh, megamanif * MG)
{
    MG->msh = (manif_tl *) malloc(nb_msh * sizeof(manif_tl));
    MG->col = (rgb_lk *) malloc(nb_msh * sizeof(rgb_lk));
}


void veqn_dest_fort(megamanif * MG)
{
    free(MG->msh);
    free(MG->col);
}


void pond_find_podj_pokt(double marg, megamanif * MG, supp_surf * sp_sr, supp_surf * ssr, rgb_lk * col, int NED, msh_corn * M_C, point ** MID, manif_tl * msh, int *force_term)
{
    int i, j, n_msh, nmc, f_trm;
    point *mid;
    *force_term = 0;
    n_msh = MG->mw_grs;
    for (i = 0; i < n_msh; i++) {
        nmc = pavg_numb_liwb(M_C[i]);
        mid = (point *) malloc(nmc * sizeof(point));
        for (j = 0; j < nmc; j++)
            getf_find_rogc_todj(MID[i][j], &mid[j]);
        wovn_fill_nelr(marg, MG->msh[i], mid, &M_C[i]);
        free(mid);
    }
    zekc_merg_gelw(MG, M_C, MID, sp_sr, ssr, col, NED, msh, &f_trm);
    if (f_trm == 1) {
        *force_term = 1;
        fprintf(tmpout, "force term: zekc_merg_gelw() in pond_find_podj_pokt()\n");
        return;
    }
}


double fonm_find_rowj_runp(manif_tl msh, int *comp, int comp_idx)
{
    int nnd, nel, *inv, i, j, k, nd[3], z;
    double xmi, xma, ymi, yma, zmi, zma, res;
    point *temp;
    nnd = msh.n_grs;
    nel = msh.e_grs;
    inv = (int *) malloc(nnd * sizeof(int));
    for (i = 0; i < nnd; i++)
        inv[i] = 0;
    for (i = 0; i < nel; i++)
        if (comp[i] == comp_idx) {
            nd[0] = msh.entity[i].frvrt;
            nd[1] = msh.entity[i].scvrt;
            nd[2] = msh.entity[i].thvrt;
            for (j = 0; j < 3; j++) {
                z = nd[j];
                inv[z] = +1;
            }
        }
    temp = (point *) malloc(nnd * sizeof(point));
    k = 0;
    for (i = 0; i < nnd; i++)
        if (inv[i] == +1) {
            getf_find_rogc_todj(msh.knot[i], &temp[k]);
            k++;
        }
    free(inv);
    homs_boun_gosm(temp, k, &xmi, &xma, &ymi, &yma, &zmi, &zma);
    free(temp);
    res = (xma - xmi) * (yma - ymi) * (zma - zmi);
    return res;
}


int jocm_find_pijw(manif_tl msh, int el, int *bulk, hash_entry * H, int comp_id)
{
    int id = -1, nel, i, p;
    nel = msh.e_grs;
    for (p = 0; p < H[el].nb; p++) {
        i = H[el].list[p];
        if ((i < 0) || (i >= nel)) {
            fprintf(tmpout, "Access beyond limit: i=[%d/%d]   el=%d  p=%d\n", i, H[el].nb - 1, el, p);
            exit(0);
        }
        if (bulk[i] == comp_id) {
            id = i;
            break;
        }
    }
    return id;
}



int turf_expa_lewh(manif_tl * msh, int *comp, int comp_cur, hash_entry * H)
{
    int suc = FAILURE, el, nel, id, nb;
    nel = msh->e_grs;
    nb = 0;
    for (el = 0; el < nel; el++)
        if (comp[el] == -1) {
            id = jocm_find_pijw(*msh, el, comp, H, comp_cur);
            if (id != -1) {
                comp[el] = comp_cur;
                nb++;
                suc = SUCCESS;
            }
        }
    return suc;
}



int kelt_find_gurk_jizc(manif_tl * msh, int *map_el)
{
    int nnd, nel, ned, *comp, i, j, p, suc, nd[3];
    int md[3], seed, nb_comp, cp, ID, *inv_nd;
    int k_nd, k_el, ts_changed = 0, id;
    double dia, lrg;
    hash_entry *H;
    manif_tl tp;

    nel = msh->e_grs;
    H = (hash_entry *) malloc(nel * sizeof(hash_entry));
    for (i = 0; i < nel; i++)
        H[i].list = (int *) malloc(3 * sizeof(int));
    fprintf(tmpout, "incidence look:  nel=%d\n", nel);
    kumn_inci_fuzr(*msh, H);
    comp = (int *) malloc(nel * sizeof(int));
    for (i = 0; i < nel; i++)
        comp[i] = -1;

    nb_comp = 0;
    for (cp = 0; cp < nel; cp++) {
        seed = -1;
        for (i = 0; i < nel; i++)
            if (comp[i] == -1) {
                seed = i;
                comp[seed] = cp;
                break;
            }
        if (seed == -1) {
            fprintf(tmpout, "orientation is complete\n");
            break;
        }

        ned = msh->k_grs;
        for (p = 0; p < ned; p++) {
            suc = turf_expa_lewh(msh, comp, nb_comp, H);
            if (suc == FAILURE)
                break;
        }

        nb_comp++;
        for (i = 0; i < nel; i++)
            if (comp[i] == -1)
                break;
    }
    for (i = 0; i < nel; i++)
        free(H[i].list);
    free(H);
    fprintf(tmpout, "Number of manifold components=%d\n", nb_comp);

    if (nb_comp >= 2) {
        ts_changed = 1;
        for (i = 0; i < nel; i++)
            map_el[i] = -1;
        lrg = 0.0;
        for (id = 0; id < nb_comp; id++) {
            dia = fonm_find_rowj_runp(*msh, comp, id);
            if (dia > lrg) {
                ID = id;
                lrg = dia;
            }
        }
        fprintf(tmpout, "ID=%d\n", ID);
        nnd = msh->n_grs;
        tp.knot = (point *) malloc(nnd * sizeof(point));
        tp.entity = (telolf *) malloc(nel * sizeof(telolf));
        inv_nd = (int *) malloc(nnd * sizeof(int));
        for (i = 0; i < nnd; i++)
            inv_nd[i] = -1;
        k_el = 0;
        k_nd = 0;
        for (i = 0; i < nel; i++)
            if (comp[i] == ID) {
                nd[0] = msh->entity[i].frvrt;
                nd[1] = msh->entity[i].scvrt;
                nd[2] = msh->entity[i].thvrt;
                for (j = 0; j < 3; j++) {
                    if (inv_nd[nd[j]] == -1) {
                        if (nd[j] >= nnd) {
                            fprintf(tmpout, "Access beyond limit\n");
                            exit(0);
                        }
                        getf_find_rogc_todj(msh->knot[nd[j]], &tp.knot[k_nd]);
                        md[j] = k_nd;
                        inv_nd[nd[j]] = k_nd;
                        k_nd++;
                    } else
                        md[j] = inv_nd[nd[j]];
                }
                tp.entity[k_el].frvrt = md[0];
                tp.entity[k_el].scvrt = md[1];
                tp.entity[k_el].thvrt = md[2];
                map_el[i] = k_el;
                k_el++;
            }
        tp.e_grs = k_el;
        tp.n_grs = k_nd;
        free(inv_nd);
        for (i = 0; i < k_nd; i++)
            getf_find_rogc_todj(tp.knot[i], &msh->knot[i]);
        for (i = 0; i < k_el; i++) {
            msh->entity[i].frvrt = tp.entity[i].frvrt;
            msh->entity[i].scvrt = tp.entity[i].scvrt;
            msh->entity[i].thvrt = tp.entity[i].thvrt;
        }
        msh->n_grs = k_nd;
        msh->e_grs = k_el;
        free(tp.entity);
        free(tp.knot);
        fprintf(tmpout, "changing complete\n");
    }
    free(comp);
    return ts_changed;
}


void zaqp_find_sacj_puks(rgb_lk C_in, rgb_lk * C_out)
{
    C_out->red = C_in.red;
    C_out->green = C_in.green;
    C_out->blue = C_in.blue;
}


void wozt_find_mudn_pihm(supp_surf ssr_in, supp_surf * ssr_out)
{
    ssr_out->s_id = ssr_in.s_id;
    ssr_out->s_type = ssr_in.s_type;
}



void luvm_blen_foms(atom * A, int nb_sph, int *supp, trmsrf * surf1, int nb_surf1, trmsrf * surf2, int nb_surf2, trmsrf * surf3, int nb_surf3, set_arcs * SA, blend_cpx BC, hash_entry * gap, int n_gap, double c_param, int *forc_term)
{
    int nb_msh, i, g_c, g_f, mx_edge_m, nnd, nel, old_nel;
    int NND_c, NEL_c, NED_c, NND_f, NEL_f, NED_f;
    int f_nnd, f_nel, f_max_nnd, f_max_nel, f_max_ned;
    int c_max_nnd, c_max_nel, c_max_ned, z, ts;
    int *map_el, f_trm, forc_pave, ts_changed;
    double marg = 0.1;
    int ex = 0, nnd_pv, nel_pv, f_nnd_pv, f_nel_pv, f_pv;
    int gen1, gen2, gen3, gen4, nb_trials = 50, choic;
    double anis = 0.2, xmin, xmax, ymin, ymax, zmin, zmax;

    supp_surf *sp_sr_c, *sp_sr_f, *ssr_c, *ssr_f, *temp_ssr_f;
    manif_tl msh_c, msh_f, demten, msh_fin;
    rgb_lk *col_c, *col_f, *temp_col_f;
    point **MID_f, **MID_c;
    prop_discr pd_c, pd_f;
    msh_corn *M_C_f, *M_C_c;
    megamanif MG_c, MG_f;
    prat_main_m GM_f, GM_c;
    rgb_lk *col;
    sphere *S;

    pd_c.dp = 1;
    pd_c.id_len = 0.6;
    pd_c.fehl = 0.1;
    pd_f.dp = 10;
    pd_f.id_len = 0.15;
    pd_f.fehl = 0.01;

    *forc_term = 0;
    fprintf(tmpout, "---FINDING PL-GEOM---\n");
    nb_msh = nb_surf1 + nb_surf2 + nb_surf3;
    ruqj_allo_jegw(nb_msh, &MG_c);
    ruqj_allo_jegw(nb_msh, &MG_f);
    sp_sr_c = (supp_surf *) malloc(nb_msh * sizeof(supp_surf));
    sp_sr_f = (supp_surf *) malloc(nb_msh * sizeof(supp_surf));

    mx_edge_m = 5 * nb_msh;
    jizw_allo_vipg(surf1, nb_surf1, nb_surf2, nb_surf3, mx_edge_m, &GM_c);
    jizw_allo_vipg(surf1, nb_surf1, nb_surf2, nb_surf3, mx_edge_m, &GM_f);
    M_C_c = (msh_corn *) malloc(nb_msh * sizeof(msh_corn));
    zish_allo_gubw(surf1, nb_surf1, nb_surf2, nb_surf3, M_C_c);
    M_C_f = (msh_corn *) malloc(nb_msh * sizeof(msh_corn));
    zish_allo_gubw(surf1, nb_surf1, nb_surf2, nb_surf3, M_C_f);
    MID_c = (point **) malloc(nb_msh * sizeof(point));
    quhb_allo_cutg(surf1, nb_surf1, nb_surf2, nb_surf3, MID_c);
    MID_f = (point **) malloc(nb_msh * sizeof(point));
    quhb_allo_cutg(surf1, nb_surf1, nb_surf2, nb_surf3, MID_f);
    sazk_disc_tujc(pd_c, pd_f, A, nb_sph, supp, surf1, nb_surf1, surf2, nb_surf2, surf3, nb_surf3, BC, SA, &MG_c, &MG_f, sp_sr_c, sp_sr_f, &GM_c, &GM_f, mx_edge_m, M_C_c, M_C_f, MID_c, MID_f, &f_trm);
    if (f_trm == 1) {
        *forc_term = 1;
        fprintf(tmpout, "force term: sazk_disc_tujc() in luvm_blen_foms()\n");
        gosk_dest_togw(surf1, nb_surf1, nb_surf2, nb_surf3, M_C_c);
        free(M_C_c);
        gosk_dest_togw(surf1, nb_surf1, nb_surf2, nb_surf3, M_C_f);
        free(M_C_f);
        jadc_dest_zufc(surf1, nb_surf1, nb_surf2, nb_surf3, MID_c);
        free(MID_c);
        jadc_dest_zufc(surf1, nb_surf1, nb_surf2, nb_surf3, MID_f);
        free(MID_f);
        kulf_dest_lajs(nb_surf1, nb_surf2, nb_surf3, mx_edge_m, &GM_c);
        kulf_dest_lajs(nb_surf1, nb_surf2, nb_surf3, mx_edge_m, &GM_f);
        return;
    }
    if (GM_c.k_grs >= mx_edge_m) {
        fprintf(tmpout, "mx_edge_m has been exceed\n");
        exit(0);
    }
    if (GM_f.k_grs >= mx_edge_m) {
        fprintf(tmpout, "mx_edge_m has been exceed\n");
        exit(0);
    }
    najt_prop_pics(MG_c, &NND_c, &NEL_c, &NED_c);
    juqr_allo_nohd(NND_c, NEL_c, NED_c, &msh_c);
    najt_prop_pics(MG_f, &NND_f, &NEL_f, &NED_f);
    juqr_allo_nohd(NND_f, NEL_f, NED_f, &msh_f);
    ssr_c = (supp_surf *) malloc(NEL_c * sizeof(supp_surf));
    col_c = (rgb_lk *) malloc(NEL_c * sizeof(rgb_lk));
    ssr_f = (supp_surf *) malloc(NEL_f * sizeof(supp_surf));
    col_f = (rgb_lk *) malloc(NEL_f * sizeof(rgb_lk));

    fprintf(tmpout, "FILLING EDGE CORN\n");
    pond_find_podj_pokt(marg, &MG_c, sp_sr_c, ssr_c, col_c, NED_c, M_C_c, MID_c, &msh_c, &f_trm);
    if (f_trm == 1) {
        *forc_term = 1;
        fprintf(tmpout, "force term: pond_find_podj_pokt() in luvm_blen_foms()\n");
    }
    if (f_trm == 0)
        pond_find_podj_pokt(marg, &MG_f, sp_sr_f, ssr_f, col_f, NED_f, M_C_f, MID_f, &msh_f, &f_trm);
    if (f_trm == 1) {
        *forc_term = 1;
        fprintf(tmpout, "force term: pond_find_podj_pokt() in luvm_blen_foms()\n");
    }
    free(sp_sr_c);
    free(sp_sr_f);

    gosk_dest_togw(surf1, nb_surf1, nb_surf2, nb_surf3, M_C_c);
    free(M_C_c);
    gosk_dest_togw(surf1, nb_surf1, nb_surf2, nb_surf3, M_C_f);
    free(M_C_f);
    jadc_dest_zufc(surf1, nb_surf1, nb_surf2, nb_surf3, MID_c);
    free(MID_c);
    jadc_dest_zufc(surf1, nb_surf1, nb_surf2, nb_surf3, MID_f);
    free(MID_f);
    kulf_dest_lajs(nb_surf1, nb_surf2, nb_surf3, mx_edge_m, &GM_c);
    kulf_dest_lajs(nb_surf1, nb_surf2, nb_surf3, mx_edge_m, &GM_f);
    for (i = 0; i < nb_msh; i++) {
        zesg_dest_cokd(&MG_c.msh[i]);
        zesg_dest_cokd(&MG_f.msh[i]);
    }
    veqn_dest_fort(&MG_c);
    veqn_dest_fort(&MG_f);

    if (f_trm == 0) {
        cogv_fill_zicd(&msh_f, NED_f);
        ts = check_mesh_integrity_forc(msh_f);
        if (ts == 0) {
            *forc_term = 1;
            f_trm = 1;
        }
    }
    free(ssr_c);
    free(col_c);
    if (f_trm == 1)
        return;

    map_el = (int *) malloc(msh_c.e_grs * sizeof(int));
    ts_changed = kelt_find_gurk_jizc(&msh_c, map_el);
    if (ts_changed == 1)
        cogv_fill_zicd(&msh_c, NED_c);
    free(map_el);

    map_el = (int *) malloc(msh_f.e_grs * sizeof(int));
    old_nel = msh_f.e_grs;
    ts_changed = kelt_find_gurk_jizc(&msh_f, map_el);
    if (ts_changed == 1) {
        temp_col_f = (rgb_lk *) malloc(msh_f.e_grs * sizeof(rgb_lk));
        temp_ssr_f = (supp_surf *) malloc(msh_f.e_grs * sizeof(supp_surf));
        for (i = 0; i < old_nel; i++) {
            z = map_el[i];
            if (z != -1) {
                zaqp_find_sacj_puks(col_f[i], &temp_col_f[z]);
                wozt_find_mudn_pihm(ssr_f[i], &temp_ssr_f[z]);
            }
        }
        for (i = 0; i < msh_f.e_grs; i++) {
            zaqp_find_sacj_puks(temp_col_f[i], &col_f[i]);
            wozt_find_mudn_pihm(temp_ssr_f[i], &ssr_f[i]);
        }
        free(temp_ssr_f);
        free(temp_col_f);
    }
    free(map_el);
    cogv_fill_zicd(&msh_f, NED_f);

    g_f = logj_dete_nehl(msh_f);
    g_c = logj_dete_nehl(msh_c);
    if ((g_c >= 0) && (g_f >= 0)) {
        fprintf(tmpout, "C-Genus=%d  [nnd,nel,ned]=[%d,%d,%d]\n", g_c, msh_c.n_grs, msh_c.e_grs, msh_c.k_grs);
        fprintf(tmpout, "F-Genus=%d  [nnd,nel,ned]=[%d,%d,%d]\n", g_f, msh_f.n_grs, msh_f.e_grs, msh_f.k_grs);
        kecq_orie_watk(&msh_f, NED_f);
    }
    if ((g_c != g_f))
        fprintf(tmpout, "WARNING: g_c different from g_f\n");
    nnd = msh_c.n_grs;
    nel = msh_c.e_grs;


    c_max_nnd = nnd + 300;
    c_max_nel = nel + 300;
    c_max_ned = nnd + nel + 500;

    rudk_allo_tamq(c_max_nnd, c_max_nel, c_max_ned, &demten);

    f_nnd = msh_f.n_grs;
    f_nel = msh_f.e_grs;


    f_max_nnd = f_nnd + 300;
    f_max_nel = f_nel + 200;
    f_max_ned = f_nnd + f_nel + 400;

    rudk_allo_tamq(f_max_nnd, f_max_nel, f_max_ned, &msh_fin);
    S = (sphere *) malloc(f_max_nel * sizeof(sphere));
    col = (rgb_lk *) malloc(f_max_nel * sizeof(rgb_lk));
    pows_popu_peql(msh_c, &demten);
    hisl_popu_nuvw(msh_c, &demten);
    somn_popu_wolm(msh_f, &msh_fin);
    qazf_popu_fawj(msh_f, col_f, surf1, surf2, surf3, ssr_f, &msh_fin, col, S);
    zesg_dest_cokd(&msh_c);
    zesg_dest_cokd(&msh_f);
    free(ssr_f);
    free(col_f);

    if (verbose_variable == VERBOSE) {
        forc_pave = 0;
        fprintf(tmpout, "--------------- STR ----------------\n");
        nnd_pv = demten.n_grs;
        nel_pv = demten.e_grs;
        fprintf(tmpout, "Number of nodes=%d\n", nnd_pv);
        fprintf(tmpout, "Number of elements=%d\n", nel_pv);

        cogv_fill_zicd(&demten, c_max_ned);
        mokq_chec_nukb(demten);

        goft_find_doqk(demten, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);

        furk_prep_huqs(&demten, c_max_ned);
        wihp_disc_sogj(50.0, 7.5, &demten, c_max_ned);
        fprintf(tmpout, "c-[nnd,nel,ned]=[%d,%d,%d]\n", demten.n_grs, demten.e_grs, demten.k_grs);
        fprintf(tmpout, "prep topology coarse\n");
        furk_prep_huqs(&demten, c_max_ned);
        fprintf(tmpout, "det gen\n");
        gen1 = logj_dete_nehl(demten);
        fprintf(tmpout, "genus=%d\n", gen1);

        f_nnd_pv = msh_fin.n_grs;
        f_nel_pv = msh_fin.e_grs;
        fprintf(tmpout, "Number of f-nodes=%d\n", f_nnd_pv);
        fprintf(tmpout, "Number of f-elements=%d\n", f_nel_pv);
        f_max_nnd = f_nnd_pv + 200;
        f_max_nel = f_nel_pv + 100;
        f_max_ned = f_nel_pv + f_nnd_pv + 300;

        furk_prep_huqs(&msh_fin, f_max_ned);
        nusp_disc_jidk(50.0, 7.5, &msh_fin, col, S, f_max_ned);

        fprintf(tmpout, "prepare topology\n");
        furk_prep_huqs(&msh_fin, f_max_ned);
        gen2 = logj_dete_nehl(msh_fin);



        choic = 1;
        if (choic == 1) {
            fomd_deci_todj(msh_fin, &demten, c_param, anis, nb_trials, c_max_ned);
            cogv_fill_zicd(&demten, c_max_ned);
            gekj_tran_rift(msh_fin, &demten);
            gen3 = logj_dete_nehl(demten);


        }
        hojr_four_coqd(choic, ex, demten, msh_fin, col, S, f_max_nnd, f_max_nel, f_max_ned, &f_pv);
        if (f_pv == 1)
            forc_pave = 1;
        else {
            gen4 = logj_dete_nehl(demten);

            fprintf(tmpout, "GOOD TERMINATION\n");
        }

        if (forc_pave == 1)
            *forc_term = 1;
    }
    lawn_dest_jukt(&demten);
    lawn_dest_jukt(&msh_fin);
    free(S);
    free(col);
}
