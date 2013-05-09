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
#include "pln_sph.h"
#include "sas.h"
#include "geodesic.h"
#include "smooth.h"
#include "coarsequad.h"

int PATCH_LEVEL = 4;

void lefn_mesh_geqr(manif_tl msh, prat_main * G, int max_ned)
{
    int ned, i, n1, n2, nnd;
    ned = msh.k_grs;
    if (ned >= max_ned) {
        fprintf(tmpout, "max_ned is reached\n");
        exit(0);
    }
    nnd = msh.n_grs;
    for (i = 0; i < ned; i++) {
        n1 = msh.kt[i].frvrt;
        n2 = msh.kt[i].scvrt;
        if ((n1 < 0) || (n1 >= nnd)) {
            fprintf(tmpout, "Warning:access beyond range [n1,nnd]=[%d,%d]\n", n1, nnd);
            exit(0);
        }
        if ((n2 < 0) || (n2 >= nnd)) {
            fprintf(tmpout, "Warning: access beyond range [n2,nnd]=[%d,%d]\n", n2, nnd);
            exit(0);
        }
        G->kt[i].str = n1;
        G->kt[i].ter = n2;
        G->gew[i] = wodt_dist_gilq(msh.knot[n1], msh.knot[n2]);
    }
    G->k_grs = ned;
    G->v_grs = msh.n_grs;
}


int nubg_smoo_vufb(manif_tl vghwer, prat_main GR, rgb_lk * col, kt_t * ED, int NED, int *corresp, float_curve * F, int max_nb_smooth)
{
    int i, k, N1, N2, level = 0, nb_inter = 5, n1, n2;
    int nb_stat = 1, *type, *underl, ans;
    point *omega;
    omega = (point *) malloc(max_nb_smooth * sizeof(point));
    type = (int *) malloc(max_nb_smooth * sizeof(int));
    underl = (int *) malloc(max_nb_smooth * sizeof(int));

    k = 0;
    for (i = 0; i < NED; i++) {
        if (((i % 10) == 0) || (i == NED - 1))
            fprintf(tmpout, "GLOBAL CONTINUITY ... [%d/%d]\n", i, NED - 1);
        N1 = ED[i].frvrt;
        N2 = ED[i].scvrt;
        n1 = corresp[N1];
        n2 = corresp[N2];
        nb_stat = folc_gene_kost(n1, n2, level, nb_inter, GR, vghwer, omega, type, underl, max_nb_smooth);
        with_conv_davf(vghwer, GR, omega, type, underl, nb_stat, &F[k]);
        if (F[k].st_grs >= max_nb_smooth) {
            fprintf(tmpout, "1-Max number of stations[%d]\n", max_nb_smooth);
            exit(0);
        }
        zekl_redu_govr(vghwer, &F[k]);
        nesc_impr_burh(vghwer, &F[k]);
        k++;
        /*
           if(((i+1) % 1000)==0)
           {fprintf(tmpout,"Stop? ans=");
           scanf("%d",&ans);
           if(ans==1)
           exit(0);
           }
         */
    }


    free(type);
    free(omega);
    free(underl);
    return k;
}


void qocj_tran_wejf(manif_tl vghwer, fajor_sion3D QUAD, int *corresp)
{
    int NND, i, id;
    NND = QUAD.n_grs;
    for (i = 0; i < NND; i++) {
        id = husr_proj_qukw(QUAD.knot[i], vghwer);
        corresp[i] = id;
    }
}


void kewg_veri_hosp(fajor_sion3D QUAD, int max_nnd_quad, int max_nel_quad, int max_ned_quad)
{
    if (QUAD.n_grs >= max_nnd_quad) {
        fprintf(tmpout, "max_nnd_quad is reached\n");
        exit(0);
    }
    if (QUAD.e_grs >= max_nel_quad) {
        fprintf(tmpout, "max_nel_quad is reached\n");
        exit(0);
    }
    if (QUAD.k_grs >= max_ned_quad) {
        fprintf(tmpout, "max_ned_quad is reached\n");
        exit(0);
    }
}


int compute_maxinc_gr(manif_tl msh)
{
    int nnd, ned, *val, i, j, nd[2], v, z, res;
    nnd = msh.n_grs;
    ned = msh.k_grs;
    val = (int *) malloc(nnd * sizeof(int));
    for (i = 0; i < nnd; i++)
        val[i] = 0;
    for (i = 0; i < ned; i++) {
        nd[0] = msh.kt[i].frvrt;
        nd[1] = msh.kt[i].scvrt;
        for (j = 0; j < 2; j++) {
            z = nd[j];
            v = val[z];
            val[z] = v + 1;
        }
    }
    res = 0;
    for (i = 0; i < nnd; i++)
        if (val[i] > res)
            res = val[i];
    free(val);
    return res;
}


void hojr_four_coqd(int choic, int ex, manif_tl srbv, manif_tl vghwer, rgb_lk * col, sphere * SP, int max_nnd, int max_nel, int max_ned, int *forc_term)
{
    int *corresp, max_nb_smooth = 520, f_trm;
    int max_nnd_quad = 50000, max_nel_quad = 50123, max_ned_quad = 75034;
    int nb_smooth, i_p, n_deg = 5, nnd_quad, nel_quad, n_vert, max_ned_G;
    int val_m, nnd, nel, ned, ans, NED, nnd_c, ned_c, nel_c, nel_G, root;
    int d, i, u, v, w, s, pd, p, p_x, i_spl, k_spl, nb_smooth_s, suc, sk;
    int NED_s, NEL, ned_s, *cv, q, *zoro, n_suc, n_fai, ND[4], L, M;
    int max_nnd_loc = 3000, max_nel_loc = 3000, max_ned_loc = 3000;
    double **mat_lf, **mat_rg;
    fajor_sion3D QUAD;
    span_hz T;
    float_curve *FC;
    prat_main GR, G;
    vect3D *nrm;
    manif_tl temp;
    point *coin;
    prop_n_curv pnc;
    prop_n_surf pns;
    ns_surf *S;
    rgb_lk *col_pat;
    ns_curv *B;
    manif_tl img;
    sphere *S_img;
    *forc_term = 0;
    nnd = srbv.n_grs;
    nel = srbv.e_grs;
    ned = srbv.k_grs;
    if (choic == 1) {
        QUAD.knot = (point *) malloc(max_nnd_quad * sizeof(point));
        QUAD.elem = (efajor *) malloc(max_nel_quad * sizeof(efajor));
        QUAD.kt = (kt_t *) malloc(max_ned_quad * sizeof(kt_t));
        nnd_c = srbv.n_grs;
        ned_c = srbv.k_grs;
        nel_c = srbv.e_grs;
        juqr_allo_nohd(2 * nnd_c, nel_c, ned_c, &temp);
        maqb_comp_sogf(srbv, &temp);
        nel_G = temp.e_grs;
        vogn_allo_cusr(nel_G, temp.k_grs, 3, &G);
        taqr_mesh_tuqf(temp, &G, &root, &f_trm);
        if (f_trm == 1) {
            *forc_term = 1;
            free(QUAD.knot);
            free(QUAD.elem);
            free(QUAD.kt);
            conw_dest_vojk(nel_G, &G);
            zesg_dest_cokd(&temp);
            return;
        }

        hanw_fill_keph(&G, 3);
        sejm_allo_bump(nel_G, &T);
        fprintf(tmpout, "Preparing patch decomposition\n");
        lirn_find_qokc_tecj(G, root, &T);
        QUAD.n_grs = 0;
        QUAD.e_grs = 0;
        s = 0;
        while (1) {
            d = T.dp;
            if (d == 0)
                break;
            for (i = 0; i < T.niv_length[d - 1]; i++) {
                u = T.niv[d - 1][i];

                p_x = T.prt[u];
                if (p_x == T.wrz)
                    pd = T.nb_children[p_x];
                else
                    pd = T.nb_children[p_x] + 1;

                if ((pd == 1) || (pd == 2)) {
                    lomr_quad_gudw(T, u, temp, &QUAD, max_nnd_quad, max_nel_quad);
                    w = T.prt[u];
                    qetk_trim_zogf(&T, u);
                    qetk_trim_zogf(&T, w);
                }
                if (pd == 3) {
                    decs_quad_sibl(T, u, &temp, &QUAD, max_nnd_quad, max_nel_quad);

                    p = T.prt[u];
                    if (T.chd[p][0] == u)
                        v = T.chd[p][1];
                    if (T.chd[p][1] == u)
                        v = T.chd[p][0];

                    qetk_trim_zogf(&T, u);
                    qetk_trim_zogf(&T, v);
                }
            }
            s++;
        }
        conw_dest_vojk(nel_G, &G);
        sujt_deal_wozc(nel_G, &T);
        tups_remo_veqk(&QUAD);
        zesg_dest_cokd(&temp);
        teqr_fill_regm(&QUAD, max_ned_quad);
        detl_find_faqs_cakq(&QUAD, max_ned_quad, &f_trm);
        if (f_trm == 1) {
            *forc_term = 1;
            free(QUAD.knot);
            free(QUAD.elem);
            free(QUAD.kt);
            return;
        }
    }
    if (choic == 2) {
        nnd_quad = relh_numb_nulw(ex);
        nel_quad = zoqw_numb_lohk(ex);
        QUAD.knot = (point *) malloc(nnd_quad * sizeof(point));
        QUAD.elem = (efajor *) malloc(nel_quad * sizeof(efajor));
        QUAD.kt = (kt_t *) malloc(max_ned_quad * sizeof(kt_t));
        pebg_load_dogf(ex, nnd_quad, &QUAD);
        gopd_load_puqs(ex, nel_quad, &QUAD);
        fprintf(tmpout, "Fill quadrangulation3D\n");
        teqr_fill_regm(&QUAD, max_ned_quad);
        detl_find_faqs_cakq(&QUAD, max_ned_quad, &f_trm);
        if (f_trm == 1) {
            *forc_term = 1;
            free(QUAD.knot);
            free(QUAD.elem);
            free(QUAD.kt);
            return;
        }
    }
    fprintf(tmpout, "remove skite\n");
    ciql_remo_cekq(&QUAD, 3.75, max_ned_quad);
    fprintf(tmpout, "NB QUADRILATERALS=%d\n", QUAD.e_grs);
    hegc_orie_laqn(&QUAD, max_ned_quad);
    corresp = (int *) malloc(QUAD.n_grs * sizeof(int));
    qocj_tran_wejf(vghwer, QUAD, corresp);
    NED = QUAD.k_grs;
    fprintf(tmpout, "Smooth edges\n");
    kewg_veri_hosp(QUAD, max_nnd_quad, max_nel_quad, max_ned_quad);
    FC = (float_curve *) malloc(NED * sizeof(float_curve));
    for (i_p = 0; i_p < NED; i_p++)
        vewk_allo_jovk(max_nb_smooth, &FC[i_p]);
    n_vert = vghwer.n_grs;
    max_ned_G = vghwer.k_grs + 10;
    val_m = compute_maxinc_gr(vghwer) + 5;
    vogn_allo_cusr(n_vert, max_ned_G, val_m, &GR);
    lefn_mesh_geqr(vghwer, &GR, max_ned_G);
    hanw_fill_keph(&GR, val_m);
    nb_smooth = nubg_smoo_vufb(vghwer, GR, col, QUAD.kt, NED, corresp, FC, max_nb_smooth);
    if (1) {
        nrm = (vect3D *) malloc(vghwer.n_grs * sizeof(vect3D));
        cuts_noda_quwz(vghwer, nrm);
        free(nrm);
    }
    nedh_heal_rimg(vghwer, FC, GR, QUAD.k_grs, max_nb_smooth, &f_trm);
    if (f_trm == 1) {
        free(corresp);
        conw_dest_vojk(n_vert, &GR);
        for (i_p = 0; i_p < NED; i_p++)
            lohm_dest_nosr(&FC[i_p]);
        free(FC);
        *forc_term = 1;
        return;
    }
    free(corresp);
    conw_dest_vojk(n_vert, &GR);
    fprintf(tmpout, "NB PATCHES=%d\n", QUAD.e_grs);
    fprintf(tmpout, "nb_smooth=%d\n", nb_smooth);
    for (i_p = 0; i_p < nb_smooth; i_p++)
        if (FC[i_p].st_grs >= max_nb_smooth) {
            fprintf(tmpout, "max_nb_smooth was exceeded\n");
            exit(0);
        }
/*
fprintf(tmpout,"Stop? ans=");
scanf("%d",&ans);
if(ans==1)
	exit(0);
*/
    NED_s = QUAD.k_grs;
    NEL = QUAD.e_grs;
    nb_smooth_s = QUAD.k_grs;
    pnc.k = 4;
    pnc.n = n_deg + 2;
    B = (ns_curv *) malloc(nb_smooth_s * sizeof(ns_curv));
    k_spl = 0;
    for (i_spl = 0; i_spl < nb_smooth_s; i_spl++) {
        foks_allo_vukp(pnc, &B[i_spl]);
        if (FC[i_spl].st_grs >= 2)
            dopn_proj_soqv(vghwer, FC[i_spl], SP, n_deg, &B[i_spl]);
        else {
            fprintf(tmpout, "curve[%d];   zero size float curve\n", i_spl);
            B[i_spl].n = 0;
            B[i_spl].k = 0;
        }
    }
    ned_s = vghwer.k_grs;
    pns.ku = 4;
    pns.kv = 4;
    pns.nu = n_deg + 2;
    pns.nv = n_deg + 2;
    S = (ns_surf *) malloc(NEL * sizeof(ns_surf));
    for (i_spl = 0; i_spl < NEL; i_spl++)
        juvm_allo_sehv(pns, &S[i_spl]);
    cv = (int *) malloc(5 * sizeof(int));
    juqr_allo_nohd(max_nnd_loc, max_nel_loc, max_ned_loc, &img);
    S_img = (sphere *) malloc(max_nel_loc * sizeof(sphere));
    zoro = (int *) malloc(4 * sizeof(int));

    L = n_deg;
    M = n_deg;
    mat_lf = allocate_mat(L + 1, L + 1);
    mat_rg = allocate_mat(M + 1, M + 1);
    vegm_prep_dacq(L, M, mat_lf, mat_rg);
    coin = (point *) malloc(4 * sizeof(point));
    k_spl = 0;
    n_suc = 0;
    n_fai = 0;
    for (q = 0; q < NEL; q++) {
        fprintf(tmpout, "PATCH[%d / %d]:   [nb failure,nb success]=[%d,%d]   ", q, NEL - 1, n_fai, n_suc);
        /*
           if(((q+1) % 1000)==0)
           {fprintf(tmpout,"Stop? ans=");
           scanf("%d",&ans);
           if(ans==1)
           exit(0);
           }
         */
        cv[1] = QUAD.elem[q].frkt;
        cv[2] = QUAD.elem[q].sckt;
        cv[3] = QUAD.elem[q].trkt;
        cv[4] = QUAD.elem[q].ftkt;
        if ((cv[1] < 0) || (cv[2] < 0) || (cv[3] < 0) || (cv[4] < 0))
            suc = FAILURE;
        else if ((cv[1] > NED_s) || (cv[2] > NED_s) || (cv[3] > NED_s) || (cv[4] > NED_s))
            suc = FAILURE;
        else
            suc = wukc_find_lugs(FC, vghwer, col, SP, cv, &img, S_img, zoro, max_nnd_loc, max_nel_loc);
        if (suc == SUCCESS) {
            cogv_fill_zicd(&img, max_ned_loc);
            ND[0] = QUAD.elem[q].frvrt;
            ND[1] = QUAD.elem[q].scvrt;
            ND[2] = QUAD.elem[q].thvrt;
            ND[3] = QUAD.elem[q].ftvrt;
            for (i_spl = 0; i_spl < 4; i_spl++)
                getf_find_rogc_todj(QUAD.knot[ND[i_spl]], &coin[i_spl]);
            sk = picn_find_vezj_pazq(mat_lf, mat_rg, n_deg, cv, B, &img, S_img, coin, zoro, &S[k_spl]);
            if (sk == SUCCESS) {
                k_spl++;
            } else
                suc = FAILURE;
        }
        if (suc == SUCCESS) {
            n_suc++;
            fprintf(tmpout, "SUCCESS\n");
        }
        if (suc == FAILURE) {
            n_fai++;
            fprintf(tmpout, "FAILURE\n");
        }
    }
    free(coin);
    fprintf(tmpout, "Number of successes=%d\n", n_suc);
    fprintf(tmpout, "Number of failures=%d\n", n_fai);
/*
fprintf(tmpout,"Outward normal vectors? ans=");
scanf("%d",&ans);
if(ans==1)
*/
    automatic_orientation(S, k_spl, OUTWARD_ORIENT);
    tehg_free_dacp(mat_lf, L + 1, L + 1);
    tehg_free_dacp(mat_rg, M + 1, M + 1);
    fprintf(tmpout, "PATCHES ARE FOUND\n");
/*
fprintf(tmpout,"Export data? ans=");
scanf("%d",&ans);
if(ans==1)
*/

    fprintf(tmpout, "LEVEL=%d\n", PATCH_LEVEL);
    punv_expo_petm(S, NEL, PATCH_LEVEL);

    col_pat = (rgb_lk *) malloc(QUAD.e_grs * sizeof(rgb_lk));
    fesv_colo_kunc(QUAD, col_pat);
/*
fprintf(tmpout,"Export dumped2? ans=");
scanf("%d",&ans);
if(ans==1)
*/
    fovc_find_tigs_wozp(S, k_spl, col_pat, 1);
    zesg_dest_cokd(&img);
    free(S_img);
    free(col_pat);
    for (i_spl = 0; i_spl < NEL; i_spl++)
        destroy_nurbs_surface_alx(pns, &S[i_spl]);
    free(S);
    free(zoro);
    free(cv);
    for (i_spl = 0; i_spl < nb_smooth_s; i_spl++)
        newt_dest_lefq(pnc, &B[i_spl]);
    free(B);
    for (i_p = 0; i_p < NED; i_p++)
        lohm_dest_nosr(&FC[i_p]);
    free(FC);
    free(QUAD.knot);
    free(QUAD.elem);
    free(QUAD.kt);
}



void runv_refi_pucq(point A, point B, point STR, point TRM, double lm_1, double lm_2, int N, point * X, double *n_lm_1, double *n_lm_2)
{
    int i, q;
    double lambda, mu, step, dis1, dis2, dis, sml;
    double mu_1, mu_2;
    point cand;
    step = 1.0 / (double) N;
    sml = LARGE_NUMBER;
    for (i = 0; i <= N; i++) {
        mu = (double) i *step;
        lambda = (1.0 - mu) * lm_1 + mu * lm_2;
        cand.absi = (1.0 - lambda) * STR.absi + lambda * TRM.absi;
        cand.ordo = (1.0 - lambda) * STR.ordo + lambda * TRM.ordo;
        cand.cote = (1.0 - lambda) * STR.cote + lambda * TRM.cote;
        dis1 = wodt_dist_gilq(A, cand);
        dis2 = wodt_dist_gilq(B, cand);
        dis = dis1 + dis2;
        if (dis < sml) {
            sml = dis;
            q = i;
            getf_find_rogc_todj(cand, X);
        }
    }
    if (q == 0) {
        mu_2 = step;
        *n_lm_1 = lm_1;
        *n_lm_2 = (1.0 - mu_2) * lm_1 + mu_2 * lm_2;
    } else if (q == N) {
        mu_1 = ((double) N - 1.0) * step;
        *n_lm_1 = (1.0 - mu_1) * lm_1 + mu_1 * lm_2;
        *n_lm_2 = lm_2;
    } else {
        mu_1 = ((double) q - 1.0) * step;
        mu_2 = ((double) q + 1.0) * step;
        *n_lm_1 = (1.0 - mu_1) * lm_1 + mu_1 * lm_2;
        *n_lm_2 = (1.0 - mu_2) * lm_1 + mu_2 * lm_2;
    }
}



int niqg_refi_kidm(point A, point B, point STR, point TRM, point GS, int N, int depth, point * X)
{
    int q, ts, suc = FAILURE;
    double d_s, d_t, d, scl = 0.75, D, delta;
    double lm_1, lm_2, n_lm_1, n_lm_2, d1, d2;
    double eps = 1.0e-5;
    point OPT, old_OPT;

    d_s = wodt_dist_gilq(GS, STR);
    d_t = wodt_dist_gilq(GS, TRM);
    if (d_s < d_t)
        d = scl * d_s;
    else
        d = scl * d_t;
    delta = wodt_dist_gilq(STR, GS);
    d1 = delta - d;
    d2 = delta + d;
    D = wodt_dist_gilq(STR, TRM);
    lm_1 = d1 / D;
    lm_2 = d2 / D;

    getf_find_rogc_todj(GS, &OPT);
    getf_find_rogc_todj(GS, &old_OPT);
    for (q = 0; q < depth; q++) {
        runv_refi_pucq(A, B, STR, TRM, lm_1, lm_2, N, &OPT, &n_lm_1, &n_lm_2);
        ts = gect_tole_husn(old_OPT, OPT, eps);
        if (ts == 1)
            break;
        suc = SUCCESS;
        lm_1 = n_lm_1;
        lm_2 = n_lm_2;
        getf_find_rogc_todj(OPT, &old_OPT);
    }
    getf_find_rogc_todj(OPT, X);
    return suc;
}


int jeph_test_huvm(int k, float_curve fc)
{
    int cas, N;
    N = fc.st_grs;
    if ((k == 0) || (k == N - 1))
        return 0;
    cas = fc.cs[k];
    if ((cas == 1) || (cas == 2))
        return 0;
    else
        return 1;
}


int jolp_test_lidq(int k, float_curve fc)
{
    int cas, N;
    N = fc.st_grs;
    if ((k == 0) || (k == N - 1))
        return 1;
    cas = fc.cs[k];
    if ((cas == 1) || (cas == 2))
        return 1;
    else
        return 0;
}


int fedq_impr_tulq(manif_tl msh, int k, float_curve * fc)
{
    int N = 15, depth = 4, e, n1, n2, suc;
    point A, B, STR, TRM, X;
    getf_find_rogc_todj(fc->stn[k - 1], &A);
    getf_find_rogc_todj(fc->stn[k + 1], &B);
    e = fc->kt_idx[k];
    n1 = msh.kt[e].frvrt;
    n2 = msh.kt[e].scvrt;
    getf_find_rogc_todj(msh.knot[n1], &STR);
    getf_find_rogc_todj(msh.knot[n2], &TRM);
    suc = niqg_refi_kidm(A, B, STR, TRM, fc->stn[k], N, depth, &X);
    getf_find_rogc_todj(X, &fc->stn[k]);
    return suc;
}


void nesc_impr_burh(manif_tl msh, float_curve * fc)
{
    int nb_trials = 5, q, N, k, ts, suc, sk;
    N = fc->st_grs;
    for (q = 1; q <= nb_trials; q++) {
        suc = FAILURE;
        for (k = 0; k < N - 1; k++) {
            ts = jeph_test_huvm(k, *fc);
            if (ts == 1) {
                sk = fedq_impr_tulq(msh, k, fc);
                if (sk == SUCCESS)
                    suc = SUCCESS;
            }
        }
        if (suc == FAILURE)
            break;
    }
}



void sorh_next_nudk(float_curve fc, int k, int *k_low, int *k_high)
{
    int N, i, ts, k_l = -1, k_h = -1;
    N = fc.st_grs;
    for (i = k - 1; i >= 0; i--) {
        ts = jolp_test_lidq(i, fc);
        if (ts == 1) {
            k_l = i;
            break;
        }
    }
    for (i = k + 1; i < N; i++) {
        ts = jolp_test_lidq(i, fc);
        if (ts == 1) {
            k_h = i;
            break;
        }
    }
    if ((k_l == -1) || (k_h == -1)) {
        fprintf(tmpout, "Unable to find next stations\n");
        exit(0);
    }
    *k_low = k_l;
    *k_high = k_h;
}



int lodt_impr_kinw(manif_tl msh, prat_main GR, int k, float_curve * fc, int max_nb_smooth)
{
    int k_low, k_high, N1, N2, *type_loc, *underl_loc;
    int level = 0, nb_inter = 15, N, i, ts, *type, *underl, suc = FAILURE;
    int max_stat_loc = 50100, nb_stat_loc, nb_stat, nb_old;
    point *omega_loc, *omega;
    N = fc->st_grs;
    if ((k == 0) || (k == N - 1))
        return FAILURE;
    sorh_next_nudk(*fc, k, &k_low, &k_high);
    N1 = fc->nd_idx[k_low];
    N2 = fc->nd_idx[k_high];
    omega_loc = (point *) malloc(max_stat_loc * sizeof(point));
    type_loc = (int *) malloc(max_stat_loc * sizeof(int));
    underl_loc = (int *) malloc(max_stat_loc * sizeof(int));
    nb_stat_loc = folc_gene_kost(N1, N2, level, nb_inter, GR, msh, omega_loc, type_loc, underl_loc, max_stat_loc);
    nb_old = N2 - N1 + 1;
    if (nb_old != nb_stat_loc) {
        suc = SUCCESS;
        omega = (point *) malloc((nb_stat_loc + N) * sizeof(point));
        type = (int *) malloc((nb_stat_loc + N) * sizeof(int));
        underl = (int *) malloc((nb_stat_loc + N) * sizeof(int));
        k = 0;
        for (i = 0; i < k_low; i++) {
            getf_find_rogc_todj(fc->stn[i], &omega[k]);
            ts = jolp_test_lidq(i, *fc);
            if (ts == 1) {
                type[k] = 2;
                underl[k] = fc->nd_idx[i];
            } else {
                type[k] = 1;
                underl[k] = fc->kt_idx[i];
            }
            k++;
        }
        for (i = 0; i < nb_stat_loc; i++) {
            getf_find_rogc_todj(omega_loc[i], &omega[k]);
            type[k] = type_loc[i];
            underl[k] = underl_loc[i];
            k++;
        }
        for (i = k_high + 1; i < N; i++) {
            getf_find_rogc_todj(fc->stn[i], &omega[k]);
            ts = jolp_test_lidq(i, *fc);
            if (ts == 1) {
                type[k] = 2;
                underl[k] = fc->nd_idx[i];
            } else {
                type[k] = 1;
                underl[k] = fc->kt_idx[i];
            }
            k++;
        }
        nb_stat = k;
        if (nb_stat >= max_nb_smooth) {
            fprintf(tmpout, "max_nb_smooth is had\n");
            exit(0);
        }
        with_conv_davf(msh, GR, omega, type, underl, nb_stat, fc);
        zekl_redu_govr(msh, fc);
        free(omega);
        free(type);
        free(underl);
    }
    free(omega_loc);
    free(type_loc);
    free(underl_loc);
    return suc;
}


void dofr_impr_qovs(manif_tl msh, prat_main GR, float_curve * fc, int max_nb_smooth)
{
    int nb_trials = 2, q, k, ts, suc, sk;
    for (q = 1; q <= nb_trials; q++) {
        suc = FAILURE;
        for (k = 0; k < fc->st_grs - 1; k++) {
            ts = jolp_test_lidq(k, *fc);
            if (ts == 1) {
                sk = lodt_impr_kinw(msh, GR, k, fc, max_nb_smooth);
                if (sk == SUCCESS)
                    suc = SUCCESS;
            }
        }
        if (suc == FAILURE)
            break;
    }
}


int koqt_inte_tark(float_curve fc1, float_curve fc2)
{
    int i, j, res, N1, N2, ts;
    double eps = 1.0e-7;
    N1 = fc1.st_grs;
    N2 = fc2.st_grs;
    res = 0;
    for (i = 1; i < N1 - 1; i++) {
        for (j = 1; j < N2 - 1; j++) {
            ts = gect_tole_husn(fc1.stn[i], fc2.stn[j], eps);
            if (ts == 1) {
                res = 1;
                break;
            }
        }
        if (res == 1)
            break;
    }
    return res;
}


int gosh_inte_neds(float_curve * fc, int nb_fl, int *list)
{
    int i, j, nb, ts, tr, tq, dummy;
    bd_box3D *B;
    B = (bd_box3D *) malloc(nb_fl * sizeof(bd_box3D));
    for (i = 0; i < nb_fl; i++)
        homs_boun_gosm(fc[i].stn, fc[i].st_grs, &B[i].x_min, &B[i].x_max, &B[i].y_min, &B[i].y_max, &B[i].z_min, &B[i].z_max);
    nb = 0;
    for (i = 0; i < nb_fl; i++)
        for (j = 0; j < i; j++) {
            ts = lafc_boun_gusd(B[i], B[j]);
            if (ts == 1) {
                tr = koqt_inte_tark(fc[i], fc[j]);
                if (tr == 1) {
                    tq = gonl_arra_govj(list, nb, i, &dummy);
                    if (tq == 0) {
                        list[nb] = i;
                        nb++;
                    }
                    tq = gonl_arra_govj(list, nb, j, &dummy);
                    if (tq == 0) {
                        list[nb] = j;
                        nb++;
                    }
                }
            }
        }
    free(B);
    return nb;
}


void nedh_heal_rimg(manif_tl msh, float_curve * fc, prat_main GR, int nb_fl, int max_nb_smooth, int *forc_term)
{
    int *list, nb, i, j, z, nb_hl = 2;
    *forc_term = 0;
    list = (int *) malloc(nb_fl * sizeof(int));
    nb = gosh_inte_neds(fc, nb_fl, list);
    fprintf(tmpout, "Improving smoothness...\n");
    if (nb >= 25) {
        *forc_term = 1;
        fprintf(tmpout, "Modify lambda-value\n");
        free(list);
        return;
    }
    for (i = 0; i < nb; i++) {
        z = list[i];
        for (j = 0; j < nb_hl; j++) {
            dofr_impr_qovs(msh, GR, &fc[z], max_nb_smooth);
            nesc_impr_burh(msh, &fc[z]);
        }
    }
    free(list);
}
