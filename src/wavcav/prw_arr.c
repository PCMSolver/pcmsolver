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
#include "eval.h"
#include "triang.h"
#include "sas.h"
#include "splinemol.h"
#include "meshsas.h"



void lurp_find_kafg_delc(double marg, int q, manif_tl * msh, msh_corn * m_c, int msh_id, point ** MID, int *forc_term)
{
    int *map_node, *map_elem, nnd, nel, cr1, cr2, pos, nnd_old, f_trm;
    msh_corn m_c_old;
    *forc_term = 0;
    pucb_allo_pakq(*m_c, &m_c_old);
    fenq_find_zoft_daqr(*m_c, &m_c_old);
    nnd = msh->n_grs;
    nel = msh->e_grs;
    map_node = (int *) malloc(nnd * sizeof(int));
    map_elem = (int *) malloc(nel * sizeof(int));
    nnd_old = msh->n_grs;

    somn_take_dogm(q, *m_c, &cr1, &cr2);

    pos = m_c->E[q].pos;
    hucz_fuse_qorw(cr1, cr2, msh, map_node, map_elem);
    munk_upda_vugz(cr1, cr2, map_node, pos, m_c);

    tuwh_fill_zupf(m_c_old, nnd_old, cr1, cr2, marg, *msh, m_c, map_node, &f_trm);
    if (f_trm == 1) {
        fprintf(tmpout, "Force term:  tuwh_fill_zupf() in lurp_find_kafg_delc()\n");
        *forc_term = 1;
    }
    free(map_node);
    free(map_elem);
    lutn_dest_mogt(&m_c_old);
}


double rusm_area_cuwl(point A, point B, point C)
{
    double res;
    vect3D T, Z, W;
    bofp_form_nukv(A, B, &T);
    bofp_form_nukv(A, C, &Z);
    cofz_cros_fits(T, Z, &W);
    res = 0.5 * biqh_norm_dapf(W);
    return res;
}


int lobj_thin_fizn(manif_tl msh, double eps)
{
    int nel, i, n1, n2, n3;
    double srf, lrg;
    nel = msh.e_grs;
    lrg = -LARGE_NUMBER;
    for (i = 0; i < nel; i++) {
        n1 = msh.entity[i].frvrt;
        n2 = msh.entity[i].scvrt;
        n3 = msh.entity[i].thvrt;
        srf = rusm_area_cuwl(msh.knot[n1], msh.knot[n2], msh.knot[n3]);
        if (srf > lrg)
            lrg = srf;
    }
    if (lrg < eps)
        return 1;
    return 0;
}



void somn_take_dogm(int q, msh_corn m_c, int *A, int *B)
{
    int a, b, ps;
    a = m_c.E[q].n1;
    b = m_c.E[q].n2;
    ps = m_c.E[q].pos;
    if (ps == -1) {
        *A = m_c.corn_ext.list[a];
        *B = m_c.corn_ext.list[b];
    } else {
        *A = m_c.corn_int[ps].list[a];
        *B = m_c.corn_int[ps].list[b];
    }
}


void wovn_fill_nelr(double marg, manif_tl msh, point * mid, msh_corn * m_c)
{
    int nin, N, nx, i, j, k, q, a[2];
    double xmi, xma, ymi, yma, zmi, zma;
    point *A;
    k = 0;

    N = m_c->corn_ext.nb;
    for (i = 0; i < N; i++) {
        nx = i + 1;
        if (nx == N)
            nx = 0;
        m_c->E[k].pos = -1;
        m_c->E[k].n1 = i;
        m_c->E[k].n2 = nx;
        getf_find_rogc_todj(mid[k], &m_c->E[k].mid_node);
        k++;
    }

    nin = m_c->h_grs;
    for (j = 0; j < nin; j++) {
        N = m_c->corn_int[j].nb;
        for (i = 0; i < N; i++) {
            nx = i + 1;
            if (nx == N)
                nx = 0;
            m_c->E[k].pos = j;
            m_c->E[k].n1 = i;
            m_c->E[k].n2 = nx;
            getf_find_rogc_todj(mid[k], &m_c->E[k].mid_node);
            k++;
        }
    }
    m_c->k_grs = k;

    A = (point *) malloc(2 * sizeof(point));
    for (q = 0; q < m_c->k_grs; q++) {
        somn_take_dogm(q, *m_c, &a[0], &a[1]);
        for (i = 0; i < 2; i++)
            getf_find_rogc_todj(msh.knot[a[i]], &A[i]);
        homs_boun_gosm(A, 2, &xmi, &xma, &ymi, &yma, &zmi, &zma);
        m_c->E[q].BB.x_min = xmi - marg;
        m_c->E[q].BB.x_max = xma + marg;

        m_c->E[q].BB.y_min = ymi - marg;
        m_c->E[q].BB.y_max = yma + marg;

        m_c->E[q].BB.z_min = zmi - marg;
        m_c->E[q].BB.z_max = zma + marg;
    }
    free(A);
}



double qehs_dist_zugt(int q1, int q2, msh_corn m_c1, msh_corn m_c2, manif_tl msh1, manif_tl msh2)
{
    int a, b;
    double dis, d1, d2;
    point A1, B1, A2, B2;
    somn_take_dogm(q1, m_c1, &a, &b);
    getf_find_rogc_todj(msh1.knot[a], &A1);
    getf_find_rogc_todj(msh1.knot[b], &B1);

    somn_take_dogm(q2, m_c2, &a, &b);
    getf_find_rogc_todj(msh2.knot[a], &A2);
    getf_find_rogc_todj(msh2.knot[b], &B2);

    d1 = wodt_dist_gilq(A1, A2) + wodt_dist_gilq(B1, B2);
    d2 = wodt_dist_gilq(A1, B2) + wodt_dist_gilq(B1, A2);
    if (d1 < d2)
        dis = d1;
    else
        dis = d2;
    return dis;
}



int tanv_asso_qupr(int p, int q, megamanif MG, int *obsol, msh_corn * M_C, int *p_asc, int *q_asc)
{
    int i, j, ts, res = FAILURE;
    double d, sml;
    sml = LARGE_NUMBER;
    for (i = 0; i < MG.mw_grs; i++)
        if ((i != p) && (obsol[i] == 0))
            for (j = 0; j < M_C[i].k_grs; j++) {
                ts = lafc_boun_gusd(M_C[p].E[q].BB, M_C[i].E[j].BB);
                if (ts == 1) {
                    d = qehs_dist_zugt(q, j, M_C[p], M_C[i], MG.msh[p], MG.msh[i]);
                    if (d < sml) {
                        sml = d;
                        *p_asc = i;
                        *q_asc = j;
                        res = SUCCESS;
                    }
                }
            }
    return res;
}


void cufq_othe_gitl(msh_corn * M_C, int p_asc, int q_asc, int *q_other)
{
    int ned, i, pos;
    *q_other = -1;
    pos = M_C[p_asc].E[q_asc].pos;
    ned = M_C[p_asc].k_grs;
    for (i = 0; i < ned; i++)
        if (i != q_asc)
            if (M_C[p_asc].E[i].pos == pos) {
                *q_other = i;
                break;
            }
}



int qevp_find_lovn_fejq(megamanif MG, msh_corn * M_C, int p_asc, int q_asc, int p, int q, int *P_ASC, int *Q_ASC, point ** MID)
{
    int pos, q_oth;
    double d1, d2;
    point REF;
    pos = M_C[p_asc].E[q_asc].pos;
    if (pos == -1) {
        if (M_C[p_asc].corn_ext.nb != 2)
            return 0;
    } else {
        if (M_C[p_asc].corn_int[pos].nb != 2)
            return 0;
    }
    cufq_othe_gitl(M_C, p_asc, q_asc, &q_oth);
    if (q_oth == -1) {
        fprintf(tmpout, "Unable to find sibling\n");
        exit(0);
    }
    getf_find_rogc_todj(M_C[p].E[q].mid_node, &REF);
    d1 = wodt_dist_gilq(REF, M_C[p_asc].E[q_asc].mid_node);
    d2 = wodt_dist_gilq(REF, M_C[p_asc].E[q_oth].mid_node);
    if (d1 < d2) {
        *P_ASC = p_asc;
        *Q_ASC = q_asc;
    } else {
        *P_ASC = p_asc;
        *Q_ASC = q_oth;
    }
    return 1;
}



int welq_asso_benv(int p, int q, megamanif MG, int *obsol, msh_corn * M_C, point ** MID, int *p_asc, int *q_asc, int **exc)
{
    int i, j, ts, res = FAILURE;
    int P_ASC, Q_ASC;
    double d, sml;
    sml = LARGE_NUMBER;
    for (i = 0; i < MG.mw_grs; i++)
        if ((i != p) && (obsol[i] == 0))
            for (j = 0; j < M_C[i].k_grs; j++)
                if (exc[i][j] == 0) {
                    ts = lafc_boun_gusd(M_C[p].E[q].BB, M_C[i].E[j].BB);
                    if (ts == 1) {
                        d = qehs_dist_zugt(q, j, M_C[p], M_C[i], MG.msh[p], MG.msh[i]);
                        if (d < sml) {
                            sml = d;
                            *p_asc = i;
                            *q_asc = j;
                            res = SUCCESS;
                        }
                    }
                }
    if (res == SUCCESS) {
        ts = qevp_find_lovn_fejq(MG, M_C, *p_asc, *q_asc, p, q, &P_ASC, &Q_ASC, MID);
        if (ts == 1) {
            fprintf(tmpout, "It has a sibling:  sml=%f\n", sml);
            *p_asc = P_ASC;
            *q_asc = Q_ASC;
        }
    }
    return res;
}



void netg_find_cofv_fuvm(conn_seq C, int a, int b, hash_entry * H, hash_entry * K)
{
    int ts_a, id_a, i, ts_b, id_b, q, N, cr, nx, k, *temp, M;
    q = -1;
    for (i = 0; i < C.nb_comp; i++) {
        ts_a = gonl_arra_govj(C.seq[i].list, C.seq[i].nb, a, &id_a);
        if (ts_a == 1) {
            ts_b = gonl_arra_govj(C.seq[i].list, C.seq[i].nb, b, &id_b);
            if (ts_b == 0) {
                fprintf(tmpout, "a and b must be on the same connected components\n");
                exit(0);
            }
            q = i;
        }
    }
    if (q == -1) {
        fprintf(tmpout, "Unable to find supporting component\n");
        exit(0);
    }
    N = C.seq[q].nb;

    cr = id_a;
    k = 0;
    for (i = 0; i < N; i++) {
        H->list[k] = C.seq[q].list[cr];
        k++;
        if (cr == id_b)
            break;
        nx = cr + 1;
        if (nx == N)
            nx = 0;
        cr = nx;
    }
    H->nb = k;

    cr = id_b;
    k = 0;
    for (i = 0; i < N; i++) {
        K->list[k] = C.seq[q].list[cr];
        k++;
        if (cr == id_a)
            break;
        nx = cr + 1;
        if (nx == N)
            nx = 0;
        cr = nx;
    }
    K->nb = k;

    M = K->nb;
    temp = (int *) malloc(M * sizeof(int));
    for (i = 0; i < M; i++)
        temp[i] = K->list[M - 1 - i];
    for (i = 0; i < M; i++)
        K->list[i] = temp[i];
    free(temp);
}



int tesq_path_hegn(manif_tl msh, int nd1, int nd2, int pos, msh_corn m_c, conn_seq C, int *pth)
{
    int ts1, ts2, id1, id2, i, q, sz, nd, ts, dummy, N;
    hash_entry H, K;
    sz = 0;
    for (i = 0; i < C.nb_comp; i++)
        sz = sz + C.seq[i].nb;
    H.list = (int *) malloc(sz * sizeof(int));
    K.list = (int *) malloc(sz * sizeof(int));
    netg_find_cofv_fuvm(C, nd1, nd2, &H, &K);

    if (pos == -1) {
        ts1 = gonl_arra_govj(m_c.corn_ext.list, m_c.corn_ext.nb, nd1, &id1);
        ts2 = gonl_arra_govj(m_c.corn_ext.list, m_c.corn_ext.nb, nd2, &id2);
        if ((ts1 == 0) || (ts2 == 0)) {
            fprintf(tmpout, "nd1 and nd2 must be on exterior boundary\n");
            exit(0);
        }
        for (i = 0; i < m_c.corn_ext.nb; i++) {
            nd = m_c.corn_ext.list[i];
            if ((nd != nd1) && (nd != nd2)) {
                q = nd;
                break;
            }
        }
    }

    if (pos >= 0) {
        ts1 = gonl_arra_govj(m_c.corn_int[pos].list, m_c.corn_int[pos].nb, nd1, &id1);
        ts2 = gonl_arra_govj(m_c.corn_int[pos].list, m_c.corn_int[pos].nb, nd2, &id2);
        if ((ts1 == 0) || (ts2 == 0)) {
            fprintf(tmpout, "nd1 and nd2 must be on exterior boundary\n");
            exit(0);
        }
        for (i = 0; i < m_c.corn_int[pos].nb; i++) {
            nd = m_c.corn_int[pos].list[i];
            if ((nd != nd1) && (nd != nd2)) {
                q = nd;
                break;
            }
        }
    }

    ts = gonl_arra_govj(H.list, H.nb, q, &dummy);
    if (ts == 0) {
        for (i = 0; i < H.nb; i++)
            pth[i] = H.list[i];
        N = H.nb;
    } else {
        for (i = 0; i < K.nb; i++)
            pth[i] = K.list[i];
        N = K.nb;
    }
    free(H.list);
    free(K.list);
    return N;
}


int socw_fami_favq(fam_nodes F, int p, int q, int *id)
{
    int res, i;
    res = 0;
    for (i = 0; i < F.n_grs; i++)
        if ((F.msh_id[i] == p) && (F.knot_id[i] == q)) {
            *id = i;
            res = 1;
            break;
        }
    return res;
}


void hadk_find_dizg_tasf(fam_nodes F, fam_nodes * G)
{
    int i, N;
    N = F.n_grs;
    getf_find_rogc_todj(F.X, &G->X);
    for (i = 0; i < N; i++) {
        G->msh_id[i] = F.msh_id[i];
        G->knot_id[i] = F.knot_id[i];
    }
    G->n_grs = N;
}


void ticb_disp_gorq(megamanif MG, fam_nodes F)
{
    int i, N, m, q;
    N = F.n_grs;
    fprintf(tmpout, "X=[%f,%f,%f]\n", F.X.absi, F.X.ordo, F.X.cote);
    for (i = 0; i < N; i++) {
        m = F.msh_id[i];
        q = F.knot_id[i];
        fprintf(tmpout, "member[%d]: msh_id=%d  node_id=%d  [%f,%f,%f]\n", i, m, q, MG.msh[m].knot[q].absi, MG.msh[m].knot[q].ordo, MG.msh[m].knot[q].cote);
    }
}


int kusd_over_wets(fam_nodes F, fam_nodes G)
{
    int dummy, ts, i, p, q, res;
    res = 0;
    for (i = 0; i < G.n_grs; i++) {
        p = G.msh_id[i];
        q = G.knot_id[i];
        ts = socw_fami_favq(F, p, q, &dummy);
        if (ts == 1) {
            res = 1;
            break;
        }
    }
    return res;
}


void nofl_merg_torg(fam_nodes F, fam_nodes G, fam_nodes * H)
{
    int dummy, ts, i, p, q, nb;
    hadk_find_dizg_tasf(F, H);
    nb = H->n_grs;
    for (i = 0; i < G.n_grs; i++) {
        p = G.msh_id[i];
        q = G.knot_id[i];
        ts = socw_fami_favq(*H, p, q, &dummy);
        if (ts == 0) {
            H->msh_id[nb] = p;
            H->knot_id[nb] = q;
            nb++;
        }
    }
    H->n_grs = nb;
    H->X.absi = 0.5 * (F.X.absi + G.X.absi);
    H->X.ordo = 0.5 * (F.X.ordo + G.X.ordo);
    H->X.cote = 0.5 * (F.X.cote + G.X.cote);
}


int gosq_unif_rukm(fam_nodes * F, int n_prw, int max_memb)
{
    int nb, ts, i, j, ind, id, N;
    fam_nodes *mrg, temp;
    mrg = (fam_nodes *) malloc(n_prw * sizeof(fam_nodes));
    for (i = 0; i < n_prw; i++) {
        mrg[i].msh_id = (int *) malloc(max_memb * sizeof(int));
        mrg[i].knot_id = (int *) malloc(max_memb * sizeof(int));
    }
    nb = 0;
    for (i = 0; i < n_prw; i++) {
        ind = 1;
        for (j = 0; j < nb; j++) {
            ts = kusd_over_wets(mrg[j], F[i]);
            if (ts == 1) {
                id = j;
                ind = 2;
                break;
            }
        }

        if (ind == 1) {
            if (F[i].n_grs >= max_memb) {
                fprintf(tmpout, "1-max_memb is reached in reunif\n");
                exit(0);
            }
            hadk_find_dizg_tasf(F[i], &mrg[nb]);
            nb++;
        }

        if (ind == 2) {
            N = mrg[id].n_grs + F[i].n_grs;
            temp.msh_id = (int *) malloc(N * sizeof(int));
            temp.knot_id = (int *) malloc(N * sizeof(int));
            nofl_merg_torg(mrg[id], F[i], &temp);
            if (N >= max_memb) {
                fprintf(tmpout, "2-max_memb is reached in reunif\n");
                exit(0);
            }
            hadk_find_dizg_tasf(temp, &mrg[id]);
            free(temp.knot_id);
            free(temp.msh_id);
        }
    }
    for (i = 0; i < nb; i++)
        hadk_find_dizg_tasf(mrg[i], &F[i]);
    for (i = 0; i < n_prw; i++) {
        free(mrg[i].msh_id);
        free(mrg[i].knot_id);
    }
    free(mrg);
    return nb;
}


int repv_find_mowl_vurw(megamanif MG, incident_gnode ** CG, int n_corn, msh_corn * M_C, fam_nodes * F, int max_fam, int max_memb)
{
    int i, j, k, p0, q0, p1, q1, n_prw, a0, b0, a1, b1;
    int temp, ts0, ts1, m, ind, N, p, q, dummy;
    double d_a, d_b, D1, D2;
    fam_nodes *f_loc;
    f_loc = (fam_nodes *) malloc(2 * sizeof(fam_nodes));
    for (i = 0; i < 2; i++) {
        f_loc[i].msh_id = (int *) malloc(2 * sizeof(int));
        f_loc[i].knot_id = (int *) malloc(2 * sizeof(int));
    }
    n_prw = 0;
    for (i = 0; i < n_corn; i++) {
        p0 = CG[i][0].msh_id;
        q0 = CG[i][0].kt_id;
        p1 = CG[i][1].msh_id;
        q1 = CG[i][1].kt_id;
        somn_take_dogm(q0, M_C[p0], &a0, &b0);
        somn_take_dogm(q1, M_C[p1], &a1, &b1);

        D1 = wodt_dist_gilq(MG.msh[p0].knot[a0], MG.msh[p1].knot[a1]);
        D2 = wodt_dist_gilq(MG.msh[p0].knot[b0], MG.msh[p1].knot[b1]);
        d_a = D1 + D2;
        D1 = wodt_dist_gilq(MG.msh[p0].knot[a0], MG.msh[p1].knot[b1]);
        D2 = wodt_dist_gilq(MG.msh[p1].knot[a1], MG.msh[p0].knot[b0]);
        d_b = D1 + D2;
        if (d_b < d_a) {
            temp = a1;
            a1 = b1;
            b1 = temp;
        }

        f_loc[0].msh_id[0] = p0;
        f_loc[0].knot_id[0] = a0;
        f_loc[0].msh_id[1] = p1;
        f_loc[0].knot_id[1] = a1;
        f_loc[0].n_grs = 2;
        f_loc[1].msh_id[0] = p0;
        f_loc[1].knot_id[0] = b0;
        f_loc[1].msh_id[1] = p1;
        f_loc[1].knot_id[1] = b1;
        f_loc[1].n_grs = 2;

        for (k = 0; k < 2; k++) {
            ind = 1;
            for (j = 0; j < n_prw; j++) {
                ts0 = socw_fami_favq(F[j], f_loc[k].msh_id[0], f_loc[k].knot_id[0], &dummy);
                ts1 = socw_fami_favq(F[j], f_loc[k].msh_id[1], f_loc[k].knot_id[1], &dummy);
                if ((ts0 == 1) && (ts1 == 0)) {
                    m = F[j].n_grs;
                    if (m >= max_memb) {
                        ticb_disp_gorq(MG, F[j]);
                        fprintf(tmpout, "max_memb is reached\n");
                        exit(0);
                    }
                    F[j].msh_id[m] = f_loc[k].msh_id[1];
                    F[j].knot_id[m] = f_loc[k].knot_id[1];
                    F[j].n_grs = m + 1;
                    ind = 2;
                    break;
                }
                if ((ts0 == 0) && (ts1 == 1)) {
                    m = F[j].n_grs;
                    if (m >= max_memb) {
                        ticb_disp_gorq(MG, F[j]);
                        fprintf(tmpout, "max_memb is reached\n");
                        exit(0);
                    }
                    F[j].msh_id[m] = f_loc[k].msh_id[0];
                    F[j].knot_id[m] = f_loc[k].knot_id[0];
                    F[j].n_grs = m + 1;
                    ind = 2;
                    break;
                }
                if ((ts0 == 1) && (ts1 == 1)) {
                    ind = 2;
                    break;
                }
            }
            if (ind == 1) {
                if (n_prw >= max_fam) {
                    fprintf(tmpout, "max_fam is reached\n");
                    exit(0);
                }
                hadk_find_dizg_tasf(f_loc[k], &F[n_prw]);
                n_prw++;
            }
        }
    }
    for (i = 0; i < 2; i++) {
        free(f_loc[i].msh_id);
        free(f_loc[i].knot_id);
    }
    free(f_loc);

    for (i = 0; i < n_prw; i++) {
        N = F[i].n_grs;
        purq_assi_sotg(0.0, 0.0, 0.0, &F[i].X);
        for (j = 0; j < N; j++) {
            p = F[i].msh_id[j];
            q = F[i].knot_id[j];
            F[i].X.absi = F[i].X.absi + MG.msh[p].knot[q].absi;
            F[i].X.ordo = F[i].X.ordo + MG.msh[p].knot[q].ordo;
            F[i].X.cote = F[i].X.cote + MG.msh[p].knot[q].cote;
        }
        F[i].X.absi = F[i].X.absi / (double) N;
        F[i].X.ordo = F[i].X.ordo / (double) N;
        F[i].X.cote = F[i].X.cote / (double) N;
    }
    return n_prw;
}


int sunc_find_puth(msh_corn m_c, int nd)
{
    int ts, dummy, i;
    ts = gonl_arra_govj(m_c.corn_ext.list, m_c.corn_ext.nb, nd, &dummy);
    if (ts == 1)
        return -1;
    for (i = 0; i < m_c.h_grs; i++) {
        ts = gonl_arra_govj(m_c.corn_int[i].list, m_c.corn_int[i].nb, nd, &dummy);
        if (ts == 1)
            return i;
    }
    return -2;
}



int pulw_find_rukq_gojn(megamanif MG, msh_corn * M_C, incident_gnode ** CG, conn_seq * C, int n_inc, fam_nodes * F, int *forc_term)
{
    int i, j, p0, q0, p1, q1, a0, b0, a1, b1, temp, sz0, sz1;
    int *pth0, pos0, *pth1, pos1, len0, len1, L, nb, n0, n1;
    double D1, D2, d_a, d_b;
    *forc_term = 0;
    fprintf(tmpout, "Curved edges\n");
    nb = 0;
    for (i = 0; i < n_inc; i++) {
        p0 = CG[i][0].msh_id;
        q0 = CG[i][0].kt_id;
        p1 = CG[i][1].msh_id;
        q1 = CG[i][1].kt_id;
        somn_take_dogm(q0, M_C[p0], &a0, &b0);
        somn_take_dogm(q1, M_C[p1], &a1, &b1);

        D1 = wodt_dist_gilq(MG.msh[p0].knot[a0], MG.msh[p1].knot[a1]);
        D2 = wodt_dist_gilq(MG.msh[p0].knot[b0], MG.msh[p1].knot[b1]);
        d_a = D1 + D2;
        D1 = wodt_dist_gilq(MG.msh[p0].knot[a0], MG.msh[p1].knot[b1]);
        D2 = wodt_dist_gilq(MG.msh[p1].knot[a1], MG.msh[p0].knot[b0]);
        d_b = D1 + D2;
        if (d_b < d_a) {
            temp = a1;
            a1 = b1;
            b1 = temp;
        }

        sz0 = 0;
        sz1 = 0;
        for (j = 0; j < C[p0].nb_comp; j++)
            sz0 = sz0 + C[p0].seq[j].nb;
        for (j = 0; j < C[p1].nb_comp; j++)
            sz1 = sz1 + C[p1].seq[j].nb;
        pth0 = (int *) malloc(sz0 * sizeof(int));
        pth1 = (int *) malloc(sz1 * sizeof(int));
        pos0 = sunc_find_puth(M_C[p0], a0);
        pos1 = sunc_find_puth(M_C[p1], a1);
        if ((pos0 == -2) || (pos1 == -2)) {
            fprintf(tmpout, "Unable to find component positions\n");
            exit(0);
        }
        len0 = tesq_path_hegn(MG.msh[p0], a0, b0, pos0, M_C[p0], C[p0], pth0);
        len1 = tesq_path_hegn(MG.msh[p1], a1, b1, pos1, M_C[p1], C[p1], pth1);
        if (len0 != len1) {
            fprintf(tmpout, "They should have the same length\n");
            *forc_term = 1;
            free(pth0);
            free(pth1);
            return 0;
        }
        L = len0 - 2;
        for (j = 0; j < L; j++) {
            F[nb].n_grs = 2;
            n0 = pth0[j + 1];
            n1 = pth1[j + 1];
            F[nb].msh_id[0] = p0;
            F[nb].knot_id[0] = n0;
            F[nb].msh_id[1] = p1;
            F[nb].knot_id[1] = n1;
            F[nb].X.absi = 0.5 * (MG.msh[p0].knot[n0].absi + MG.msh[p1].knot[n1].absi);
            F[nb].X.ordo = 0.5 * (MG.msh[p0].knot[n0].ordo + MG.msh[p1].knot[n1].ordo);
            F[nb].X.cote = 0.5 * (MG.msh[p0].knot[n0].cote + MG.msh[p1].knot[n1].cote);
            nb++;
        }
        free(pth0);
        free(pth1);
    }
    return nb;
}


int dajh_wors_juvs(manif_tl msh)
{
    int i, ned, nb;
    ned = msh.k_grs;
    nb = 0;
    for (i = 0; i < ned; i++)
        if (msh.kt[i].scent == -1)
            nb++;
    return nb;
}



void qujs_loca_wilr(fam_nodes * F_corn, int n_prw, fam_nodes * F_edge, int n_edge, megamanif MG, int pt, supp_surf * sp_sr, supp_surf * ssr, manif_tl * msh)
{
    int nnd, ned, *ind, i, j, n1, n2, n3, p, q, *map_node, M;
    nnd = MG.msh[pt].n_grs;
    ind = (int *) malloc(nnd * sizeof(int));
    for (i = 0; i < nnd; i++)
        ind[i] = -1;
    ned = MG.msh[pt].k_grs;
    for (i = 0; i < ned; i++)
        if (MG.msh[pt].kt[i].scent == -1) {
            n1 = MG.msh[pt].kt[i].frvrt;
            n2 = MG.msh[pt].kt[i].scvrt;
            ind[n1] = +1;
            ind[n2] = +1;
        }

    map_node = (int *) malloc(nnd * sizeof(int));
    p = msh->n_grs;
    for (i = 0; i < nnd; i++)
        if (ind[i] == -1) {
            getf_find_rogc_todj(MG.msh[pt].knot[i], &msh->knot[p]);
            map_node[i] = p;
            p++;
        }
    free(ind);

    for (i = 0; i < n_prw; i++) {
        M = F_corn[i].n_grs;
        for (j = 0; j < M; j++)
            if (F_corn[i].msh_id[j] == pt) {
                q = F_corn[i].knot_id[j];
                map_node[q] = i;
            }
    }
    for (i = 0; i < n_edge; i++) {
        M = F_edge[i].n_grs;
        for (j = 0; j < M; j++)
            if (F_edge[i].msh_id[j] == pt) {
                q = F_edge[i].knot_id[j];
                map_node[q] = i + n_prw;
            }
    }
    q = msh->e_grs;
    for (i = 0; i < MG.msh[pt].e_grs; i++) {
        n1 = MG.msh[pt].entity[i].frvrt;
        n2 = MG.msh[pt].entity[i].scvrt;
        n3 = MG.msh[pt].entity[i].thvrt;
        msh->entity[q].frvrt = map_node[n1];
        msh->entity[q].scvrt = map_node[n2];
        msh->entity[q].thvrt = map_node[n3];
        ssr[q].s_id = sp_sr[pt].s_id;
        ssr[q].s_type = sp_sr[pt].s_type;
        q++;
    }
    free(map_node);
    msh->n_grs = p;
    msh->e_grs = q;
}


void sukt_comb_bejq(megamanif MG, int *obsol, fam_nodes * F_corn, int n_prw, fam_nodes * F_edge, int n_edge, supp_surf * sp_sr, supp_surf * ssr, manif_tl * msh)
{
    int i, k, p, n_msh;
    k = 0;
    for (i = 0; i < n_prw; i++) {
        getf_find_rogc_todj(F_corn[i].X, &msh->knot[k]);
        k++;
    }
    for (i = 0; i < n_edge; i++) {
        getf_find_rogc_todj(F_edge[i].X, &msh->knot[k]);
        k++;
    }
    msh->n_grs = k;
    msh->e_grs = 0;
    n_msh = MG.mw_grs;
    for (p = 0; p < n_msh; p++)
        if (obsol[p] == 0)
            qujs_loca_wilr(F_corn, n_prw, F_edge, n_edge, MG, p, sp_sr, ssr, msh);
}


void zekc_merg_gelw(megamanif * MG, msh_corn * M_C, point ** MID, supp_surf * sp_sr, supp_surf * ssr, rgb_lk * col, int max_ned, manif_tl * msh, int *forc_term)
{
    int *obsol, n_msh, i, j, k, n_inc, max_conn = 10, nb_w, nb;
    int max_fam, max_memb = 20, n_prw, n_edge, N, new_prw, f_trm;
    double eps_area = 1.0e-5, eps_fuse = 1.0e-4, rd, gr, bl, marg = 0.1;
    conn_seq *C;
    fam_nodes *F_corn, *F_edge;
    incident_gnode **CG;
    *forc_term = 0;
    n_msh = MG->mw_grs;
    obsol = (int *) malloc(n_msh * sizeof(int));
    for (i = 0; i < n_msh; i++)
        obsol[i] = 0;
    tuwq_disc_nafp(marg, M_C, MID, MG, obsol, eps_area, eps_fuse, &f_trm);
    if (f_trm == 1) {
        *forc_term = 1;
        fprintf(tmpout, "force term: tuwq_disc_nafp() in zekc_merg_gelw()\n");
        free(obsol);
        return;
    }
    for (i = 0; i < n_msh; i++)
        capn_fill_fond(&MG->msh[i]);
    C = (conn_seq *) malloc(n_msh * sizeof(conn_seq));
    fprintf(tmpout, "Find connected components\n");
    for (i = 0; i < n_msh; i++)
        if (obsol[i] == 0) {
            nb_w = dajh_wors_juvs(MG->msh[i]);
            C[i].seq = (hash_entry *) malloc(max_conn * sizeof(hash_entry));
            for (j = 0; j < max_conn; j++)
                C[i].seq[j].list = (int *) malloc(nb_w * sizeof(int));
            hojp_boun_hapw(MG->msh[i], &C[i], max_conn);
        }
    nb = 0;
    for (i = 0; i < n_msh; i++)
        nb = nb + M_C[i].k_grs;
    CG = (incident_gnode **) malloc(nb * sizeof(incident_gnode));
    for (i = 0; i < nb; i++)
        CG[i] = (incident_gnode *) malloc(2 * sizeof(incident_gnode));
    fprintf(tmpout, "Finding the incidence\n");
    n_inc = wuvp_inci_tejg(*MG, obsol, M_C, MID, CG);
    max_fam = 10 * n_msh;
    F_corn = (fam_nodes *) malloc(max_fam * sizeof(fam_nodes));
    for (i = 0; i < max_fam; i++) {
        F_corn[i].msh_id = (int *) malloc(max_memb * sizeof(int));
        F_corn[i].knot_id = (int *) malloc(max_memb * sizeof(int));
    }
    fprintf(tmpout, "Pairwise different nodes\n");
    n_prw = repv_find_mowl_vurw(*MG, CG, n_inc, M_C, F_corn, max_fam, max_memb);
    new_prw = gosq_unif_rukm(F_corn, n_prw, max_memb);
    fprintf(tmpout, "n_prw=%d   new_prw=%d\n", n_prw, new_prw);
    n_prw = new_prw;
    N = 0;
    for (i = 0; i < n_msh; i++)
        if (obsol[i] == 0)
            for (j = 0; j < C[i].nb_comp; j++)
                N = N + C[i].seq[j].nb;
    fprintf(tmpout, "prw diff curved edges\n");
    F_edge = (fam_nodes *) malloc(N * sizeof(fam_nodes));
    for (i = 0; i < N; i++) {
        F_edge[i].msh_id = (int *) malloc(2 * sizeof(int));
        F_edge[i].knot_id = (int *) malloc(2 * sizeof(int));
    }
    n_edge = pulw_find_rukq_gojn(*MG, M_C, CG, C, n_inc, F_edge, &f_trm);
    if (f_trm == 1) {
        *forc_term = 1;
        fprintf(tmpout, "force term: pulw_find_rukq_gojn() in zekc_merg_gelw()\n");
        for (i = 0; i < N; i++) {
            free(F_edge[i].msh_id);
            free(F_edge[i].knot_id);
        }
        free(F_edge);
        for (i = 0; i < nb; i++)
            free(CG[i]);
        free(CG);
        for (i = 0; i < n_msh; i++)
            if (obsol[i] == 0) {
                for (j = 0; j < max_conn; j++)
                    free(C[i].seq[j].list);
                free(C[i].seq);
            }
        free(C);
        for (i = 0; i < max_fam; i++) {
            free(F_corn[i].msh_id);
            free(F_corn[i].knot_id);
        }
        free(F_corn);
        free(obsol);
        return;
    }
    fprintf(tmpout, "n_edge=%d\n", n_edge);
    for (i = 0; i < nb; i++)
        free(CG[i]);
    free(CG);
    for (i = 0; i < n_msh; i++)
        if (obsol[i] == 0) {
            for (j = 0; j < max_conn; j++)
                free(C[i].seq[j].list);
            free(C[i].seq);
        }
    free(C);
    fprintf(tmpout, "Combination\n");
    sukt_comb_bejq(*MG, obsol, F_corn, n_prw, F_edge, n_edge, sp_sr, ssr, msh);
    fprintf(tmpout, "Assign colors\n");
    k = 0;
    for (i = 0; i < MG->mw_grs; i++)
        if (obsol[i] == 0) {
            tesr_colo_donr(i, &rd, &gr, &bl);
            for (j = 0; j < MG->msh[i].e_grs; j++) {
                col[k].red = rd;
                col[k].green = gr;
                col[k].blue = bl;
                k++;
            }
        }
    for (i = 0; i < N; i++) {
        free(F_edge[i].msh_id);
        free(F_edge[i].knot_id);
    }
    free(F_edge);
    for (i = 0; i < max_fam; i++) {
        free(F_corn[i].msh_id);
        free(F_corn[i].knot_id);
    }
    free(F_corn);
    free(obsol);
    cogv_fill_zicd(msh, max_ned);
    for (i = 0; i < msh->k_grs; i++)
        if (msh->kt[i].scent == -1) {
            fprintf(tmpout, "EL2=-1\n");
            exit(0);
        }
    fprintf(tmpout, "CLOSED\n");
    fprintf(tmpout, "COMPLETION\n");
}
