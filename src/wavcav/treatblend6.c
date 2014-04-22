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


int rusl_find_jons(manif_tl msh, int nd, int *mu, int nb, int *exc)
{
    int i, e, E, n1, n2;
    for (i = 0; i < nb; i++)
        if (exc[i] == 0) {
            E = mu[i];
            n1 = msh.kt[E].frvrt;
            n2 = msh.kt[E].scvrt;
            if ((n1 == nd) || (n2 == nd)) {
                e = E;
                exc[i] = 1;
                break;
            }
        }
    return e;
}



int jerl_conn_bukc(manif_tl msh, int *mu, int nb, int seed, int *seq, int *exc)
{
    int i, e, n1, n2, nx, N, beg, ter;
    seq[0] = seed;
    N = 1;
    for (i = 1; i < nb; i++) {
        e = rusl_find_jons(msh, seq[i - 1], mu, nb, exc);
        n1 = msh.kt[e].frvrt;
        n2 = msh.kt[e].scvrt;
        if (n1 == seq[i - 1])
            nx = n2;
        if (n2 == seq[i - 1])
            nx = n1;
        if (nx == seed)
            break;
        seq[N] = nx;
        N++;
    }
    beg = seq[0];
    ter = seq[N - 1];
    for (i = 0; i < nb; i++)
        if (exc[i] == 0) {
            e = mu[i];
            n1 = msh.kt[e].frvrt;
            n2 = msh.kt[e].scvrt;
            if (((n1 == beg) && (n2 == ter)) || ((n2 == beg) && (n1 == ter)))
                exc[i] = 1;
        }
    return N;
}



void hojp_boun_hapw(manif_tl msh, conn_seq * C, int max_conn)
{
    int *mu, nb, ned, i, *exc, p, seed, n_comp, e;
    ned = msh.k_grs;
    mu = (int *) malloc(ned * sizeof(int));
    nb = 0;
    for (i = 0; i < ned; i++)
        if (msh.kt[i].scent == -1) {
            mu[nb] = i;
            nb++;
        }
    exc = (int *) malloc(nb * sizeof(int));
    for (i = 0; i < nb; i++)
        exc[i] = 0;
    n_comp = 0;
    for (p = 0; p < nb; p++) {
        seed = -1;
        for (i = 0; i < nb; i++)
            if (exc[i] == 0) {
                e = mu[i];
                seed = msh.kt[e].frvrt;
                break;
            }
        if (seed == -1)
            break;
        if (n_comp >= max_conn) {
            fprintf(tmpout, "max_conn is reached\n");
            exit(0);
        }
        C->seq[p].nb = jerl_conn_bukc(msh, mu, nb, seed, C->seq[p].list, exc);
        n_comp++;
    }
    free(exc);
    free(mu);
    C->nb_comp = n_comp;
}


void hefn_veri_sezm(megamanif MG, msh_corn * M_C, int p, int q, int p_asc, int q_asc)
{
    int a, b, a_asc, b_asc;
    double d;
    somn_take_dogm(q, M_C[p], &a, &b);
    somn_take_dogm(q_asc, M_C[p_asc], &a_asc, &b_asc);
    d = qehs_dist_zugt(q, q_asc, M_C[p], M_C[p_asc], MG.msh[p], MG.msh[p_asc]);
    if (d > 0.05) {
        fprintf(tmpout, "a=[%f,%f,%f]\n", MG.msh[p].knot[a].absi, MG.msh[p].knot[a].ordo, MG.msh[p].knot[a].cote);
        fprintf(tmpout, "b=[%f,%f,%f]\n", MG.msh[p].knot[b].absi, MG.msh[p].knot[b].ordo, MG.msh[p].knot[b].cote);
        fprintf(tmpout, "a_asc=[%f,%f,%f]\n", MG.msh[p_asc].knot[a_asc].absi, MG.msh[p_asc].knot[a_asc].ordo, MG.msh[p_asc].knot[a_asc].cote);
        fprintf(tmpout, "b_asc=[%f,%f,%f]\n", MG.msh[p_asc].knot[b_asc].absi, MG.msh[p_asc].knot[b_asc].ordo, MG.msh[p_asc].knot[b_asc].cote);
        fprintf(tmpout, "Gap=%f\n", d);
        exit(0);
    }
}


void fagz_veri_fojm(megamanif MG, msh_corn * M_C, incident_gnode ** CG, int n_inc)
{
    int i, p, q, p_asc, q_asc;
    for (i = 0; i < n_inc; i++) {
        p = CG[i][0].msh_id;
        q = CG[i][0].kt_id;
        p_asc = CG[i][1].msh_id;
        q_asc = CG[i][1].kt_id;
        hefn_veri_sezm(MG, M_C, p, q, p_asc, q_asc);
    }
    fprintf(tmpout, "GOOD INC GNODE\n");
}



int wuvp_inci_tejg(megamanif MG, int *obsol, msh_corn * m_c, point ** MID, incident_gnode ** CG)
{
    int **exc, n_msh, i, j, ned, suc, p_asc, q_asc, p, q, nb;
    n_msh = MG.mw_grs;
    exc = (int **) malloc(n_msh * sizeof(int *));
    for (i = 0; i < n_msh; i++) {
        ned = m_c[i].k_grs;
        exc[i] = (int *) malloc(ned * sizeof(int));
        for (j = 0; j < ned; j++)
            exc[i][j] = 0;
    }
    nb = 0;
    for (p = 0; p < n_msh; p++)
        if (obsol[p] == 0) {
            ned = m_c[p].k_grs;
            for (q = 0; q < ned; q++)
                if (exc[p][q] == 0) {
                    suc = welq_asso_benv(p, q, MG, obsol, m_c, MID, &p_asc, &q_asc, exc);
                    if (suc == FAILURE) {
                        fprintf(tmpout, "No assoc\n");
                        exit(0);
                    }

                    exc[p][q] = 1;
                    exc[p_asc][q_asc] = 1;
                    CG[nb][0].msh_id = p;
                    CG[nb][0].kt_id = q;
                    CG[nb][1].msh_id = p_asc;
                    CG[nb][1].kt_id = q_asc;
                    nb++;
                }
        }
    for (i = 0; i < n_msh; i++)
        free(exc[i]);
    free(exc);

    return nb;
}



int tamb_tigh_joqm(msh_corn m_c, manif_tl msh, double eps, int *q1, int *q2)
{
    int q, a, b, ts, *st, res;
    point A, B;
    st = (int *) malloc(4 * sizeof(int));
    for (q = 0; q < 4; q++) {
        somn_take_dogm(q, m_c, &a, &b);
        getf_find_rogc_todj(msh.knot[a], &A);
        getf_find_rogc_todj(msh.knot[b], &B);
        ts = gect_tole_husn(A, B, eps);
        st[q] = 0;
        if (ts == 1)
            st[q] = 1;
    }
    res = 0;
    if ((st[0] == 1) && (st[2] == 1)) {
        res = 1;
        *q1 = 0;
        *q2 = 2;
    } else if ((st[1] == 1) && (st[3] == 1)) {
        res = 1;
        *q1 = 1;
        *q2 = 3;
    }
    free(st);
    return res;
}



void tuwq_disc_nafp(double marg, msh_corn * M_C, point ** MID, megamanif * MG, int *obsol, double eps_area, double eps_fuse, int *forc_term)
{
    int p, j, q[2], ts, tr, suc, p_res, q_res, f_trm;
    *forc_term = 0;
    for (p = 0; p < MG->mw_grs; p++)
        if ((M_C[p].h_grs == 0) && (M_C[p].corn_ext.nb == 4)) {
            ts = lobj_thin_fizn(MG->msh[p], eps_area);
            if (ts == 1) {
                tr = tamb_tigh_joqm(M_C[p], MG->msh[p], eps_fuse, &q[0], &q[1]);
                if (tr == 1) {
                    obsol[p] = 1;
                    for (j = 0; j < 2; j++) {
                        suc = tanv_asso_qupr(p, q[j], *MG, obsol, M_C, &p_res, &q_res);
                        if (suc == FAILURE) {
                            fprintf(tmpout, "Unable to find assoc\n");
                            exit(0);
                        }
                        lurp_find_kafg_delc(marg, q_res, &MG->msh[p_res], &M_C[p_res], p_res, MID, &f_trm);
                        if (f_trm == 1) {
                            *forc_term = 1;
                            fprintf(tmpout, "force term: lurp_find_kafg_delc() in tuwq_disc_nafp()\n");
                            return;
                        }
                    }
                }
            }
        }
    for (j = 0; j < MG->mw_grs; j++)
        if (obsol[j] == 1)
            fprintf(tmpout, "obsolete[%d/%d]\n", j, MG->mw_grs - 1);
}



void hucz_fuse_qorw(int cr1, int cr2, manif_tl * msh, int *map_node, int *map_elem)
{
    int nnd, nel, i, k, n1, n2, n3, im1, im2, im3;
    point X, *temp;
    X.absi = 0.5 * (msh->knot[cr1].absi + msh->knot[cr2].absi);
    X.ordo = 0.5 * (msh->knot[cr1].ordo + msh->knot[cr2].ordo);
    X.cote = 0.5 * (msh->knot[cr1].cote + msh->knot[cr2].cote);
    nnd = msh->n_grs;
    temp = (point *) malloc((nnd - 1) * sizeof(point));
    k = 0;
    for (i = 0; i < nnd; i++)
        if ((i != cr1) && (i != cr2)) {
            getf_find_rogc_todj(msh->knot[i], &temp[k]);
            map_node[i] = k;
            k++;
        }
    getf_find_rogc_todj(X, &temp[k]);
    map_node[cr1] = k;
    map_node[cr2] = k;
    k++;
    for (i = 0; i < nnd - 1; i++)
        getf_find_rogc_todj(temp[i], &msh->knot[i]);
    free(temp);
    msh->n_grs = nnd - 1;

    nel = msh->e_grs;
    k = 0;
    for (i = 0; i < nel; i++) {
        n1 = msh->entity[i].frvrt;
        n2 = msh->entity[i].scvrt;
        n3 = msh->entity[i].thvrt;
        im1 = map_node[n1];
        im2 = map_node[n2];
        im3 = map_node[n3];
        if ((im1 != im2) && (im1 != im3) && (im2 != im3)) {
            msh->entity[k].frvrt = im1;
            msh->entity[k].scvrt = im2;
            msh->entity[k].thvrt = im3;
            map_elem[i] = k;
            k++;
        } else
            map_elem[i] = -1;
    }
    msh->e_grs = k;
}



void munk_upda_vugz(int cr1, int cr2, int *map_node, int pos, msh_corn * M_C)
{
    int *tp, N, i, k, z, nin, qos, temp;
    if (pos == -1) {
        N = M_C->corn_ext.nb;
        tp = (int *) malloc(N * sizeof(int));
        k = 0;
        for (i = 0; i < N; i++)
            if (M_C->corn_ext.list[i] != cr2) {
                z = M_C->corn_ext.list[i];
                tp[k] = map_node[z];
                k++;
            }
        for (i = 0; i < k; i++)
            M_C->corn_ext.list[i] = tp[i];
        M_C->corn_ext.nb = k;
        free(tp);
    }
    if (pos >= 0) {
        N = M_C->corn_int[pos].nb;
        tp = (int *) malloc(N * sizeof(int));
        k = 0;
        for (i = 0; i < N; i++)
            if (M_C->corn_int[pos].list[i] != cr2) {
                z = M_C->corn_int[pos].list[i];
                tp[k] = map_node[z];
                k++;
            }
        for (i = 0; i < k; i++)
            M_C->corn_int[pos].list[i] = tp[i];
        M_C->corn_int[pos].nb = k;
        free(tp);
    }

    nin = M_C->h_grs;
    for (qos = pos + 1; qos < nin; qos++) {
        for (i = 0; i < M_C->corn_int[qos].nb; i++) {
            temp = M_C->corn_int[qos].list[i];
            M_C->corn_int[qos].list[i] = temp - 2;
        }
    }

}


void qonc_find_vofd_ralp(hash_entry E, hash_entry * F)
{
    int i, nb;
    nb = E.nb;
    for (i = 0; i < nb; i++)
        F->list[i] = E.list[i];
    F->nb = nb;
}


void kumh_find_cawj_fovl(kt_corn E, kt_corn * F)
{
    F->pos = E.pos;
    F->n1 = E.n1;
    F->n2 = E.n2;
    guwv_find_dagt_hujw(E.BB, &F->BB);
    getf_find_rogc_todj(E.mid_node, &F->mid_node);
}


void fenq_find_zoft_daqr(msh_corn m_c, msh_corn * M_C)
{
    int i, nb_h, nb_e;
    nb_h = m_c.h_grs;
    nb_e = m_c.k_grs;
    qonc_find_vofd_ralp(m_c.corn_ext, &M_C->corn_ext);
    for (i = 0; i < nb_h; i++)
        qonc_find_vofd_ralp(m_c.corn_int[i], &M_C->corn_int[i]);
    for (i = 0; i < nb_e; i++)
        kumh_find_cawj_fovl(m_c.E[i], &M_C->E[i]);
    M_C->h_grs = nb_h;
    M_C->k_grs = nb_e;
}


int dejh_find_rohf(int p1, int p2, int n1_old, int n2_old, int cr1, int cr2)
{
    int ind = 2;
    if ((p1 == n1_old) && (p2 == n2_old))
        ind = 1;
    if ((p2 == n1_old) && (p1 == n2_old))
        ind = 1;
    if ((cr1 == n1_old) && (p2 == n2_old))
        ind = 1;
    if ((cr2 == n1_old) && (p2 == n2_old))
        ind = 1;
    if ((p1 == n1_old) && (cr1 == n2_old))
        ind = 1;
    if ((p1 == n1_old) && (cr2 == n2_old))
        ind = 1;
    return ind;
}



int rukl_find_sewj_hodc(int cr1, int cr2, msh_corn M_C, int *inv_node, int n1_new, int n2_new)
{
    int e = -1, n1_old, n2_old, i, ned, m1, m2, pos, p1, p2, ind;
    n1_old = inv_node[n1_new];
    n2_old = inv_node[n2_new];
    ned = M_C.k_grs;
    for (i = 0; i < ned; i++) {
        pos = M_C.E[i].pos;
        m1 = M_C.E[i].n1;
        m2 = M_C.E[i].n2;
        if (pos == -1) {
            p1 = M_C.corn_ext.list[m1];
            p2 = M_C.corn_ext.list[m2];
            ind = dejh_find_rohf(p1, p2, n1_old, n2_old, cr1, cr2);
            if (ind == 1) {
                e = i;
                break;
            }
        }

        if (pos >= 0) {
            p1 = M_C.corn_int[pos].list[m1];
            p2 = M_C.corn_int[pos].list[m2];
            ind = dejh_find_rohf(p1, p2, n1_old, n2_old, cr1, cr2);
            if (ind == 1) {
                e = i;
                break;
            }
        }
    }
    return e;
}


void hift_disp_tevr(msh_corn m_c)
{
    int i, j;
    fprintf(tmpout, "nb holes=%d\n", m_c.h_grs);
    fprintf(tmpout, "nb edges=%d\n", m_c.k_grs);
    fprintf(tmpout, "Exterior corners:\n");
    for (i = 0; i < m_c.corn_ext.nb; i++)
        fprintf(tmpout, "\tcorn_ext[%d]=%d\n", i, m_c.corn_ext.list[i]);
    fprintf(tmpout, "Interior corners:\n");
    for (j = 0; j < m_c.h_grs; j++) {
        for (i = 0; i < m_c.corn_int[j].nb; i++)
            fprintf(tmpout, "\tcorn_int[%d][%d]=%d\n", j, i, m_c.corn_int[j].list[i]);
    }
    for (i = 0; i < m_c.k_grs; i++)
        fprintf(tmpout, "edge[%d]=[%d,%d]  pos=%d\n", i, m_c.E[i].n1, m_c.E[i].n2, m_c.E[i].pos);
}


void pucb_allo_pakq(msh_corn m_c, msh_corn * m_c_old)
{
    int nh, ne, i;
    nh = m_c.h_grs;
    ne = m_c.k_grs;
    m_c_old->corn_ext.list = (int *) malloc(m_c.corn_ext.nb * sizeof(int));
    m_c_old->corn_int = (hash_entry *) malloc(nh * sizeof(hash_entry));
    for (i = 0; i < nh; i++)
        m_c_old->corn_int[i].list = (int *) malloc(m_c.corn_int[i].nb * sizeof(int));
    m_c_old->E = (kt_corn *) malloc(ne * sizeof(kt_corn));
    m_c_old->h_grs = nh;
}


void lutn_dest_mogt(msh_corn * m_c_old)
{
    int nh, i;
    nh = m_c_old->h_grs;
    free(m_c_old->corn_ext.list);
    for (i = 0; i < nh; i++)
        free(m_c_old->corn_int[i].list);
    free(m_c_old->corn_int);
    free(m_c_old->E);
}



void tuwh_fill_zupf(msh_corn m_c_old, int nnd_old, int cr1, int cr2, double marg, manif_tl msh, msh_corn * m_c, int *map_node, int *proc_term)
{
    int nin, N, nx, i, j, k, q, a[2], e;
    int n1_new, n2_new, *inv_node, s;
    double xmi, xma, ymi, yma, zmi, zma;
    point *A;

    *proc_term = 0;
    inv_node = (int *) malloc(nnd_old * sizeof(int));
    for (i = 0; i < nnd_old; i++) {
        s = map_node[i];
        inv_node[s] = i;
    }
    k = 0;

    N = m_c->corn_ext.nb;
    for (i = 0; i < N; i++) {
        nx = i + 1;
        if (nx == N)
            nx = 0;
        m_c->E[k].pos = -1;
        m_c->E[k].n1 = i;
        m_c->E[k].n2 = nx;
        n1_new = m_c->corn_ext.list[i];
        n2_new = m_c->corn_ext.list[nx];
        e = rukl_find_sewj_hodc(cr1, cr2, m_c_old, inv_node, n1_new, n2_new);
        if (e == -1) {
            fprintf(tmpout, "1-Unable to find old edge\n");
            *proc_term = 1;
            free(inv_node);
            return;
        }
        getf_find_rogc_todj(m_c_old.E[e].mid_node, &m_c->E[k].mid_node);
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
            n1_new = m_c->corn_int[j].list[i];
            n2_new = m_c->corn_int[j].list[nx];
            e = rukl_find_sewj_hodc(cr1, cr2, m_c_old, inv_node, n1_new, n2_new);
            if (e == -1) {
                fprintf(tmpout, "[n1_new,n2_new]=[%d,%d]\n", n1_new, n2_new);
                fprintf(tmpout, "2-Unable to find old edge\n");
                *proc_term = 1;
                free(inv_node);
                return;
            }
            getf_find_rogc_todj(m_c_old.E[e].mid_node, &m_c->E[k].mid_node);
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
    free(inv_node);
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

