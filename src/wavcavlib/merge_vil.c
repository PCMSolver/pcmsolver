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
#include <stdlib.h>
#include "cavity.h"
#include "splinemol.h"
#include "geodesic.h"
#include "pln_sph.h"
#include "smooth.h"


void vunm_comm_dilp(fajor_sion3D quad, int E, int F, int *q1, int *q2)
{
    int ed[4], fd[4], i, k, cm[4], ts, dummy;
    ed[0] = quad.elem[E].frkt;
    ed[1] = quad.elem[E].sckt;
    ed[2] = quad.elem[E].trkt;
    ed[3] = quad.elem[E].ftkt;
    fd[0] = quad.elem[F].frkt;
    fd[1] = quad.elem[F].sckt;
    fd[2] = quad.elem[F].trkt;
    fd[3] = quad.elem[F].ftkt;
    k = 0;
    for (i = 0; i < 4; i++) {
        ts = gonl_arra_govj(fd, 4, ed[i], &dummy);
        if (ts == 1) {
            cm[k] = ed[i];
            k++;
        }
    }
    if (k != 2) {
        fprintf(tmpout, "E and F must share two edges\n");
        exit(0);
    }
    *q1 = cm[0];
    *q2 = cm[1];
}


void dezq_surr_qong(fajor_sion3D quad, int E, int F, int q1, int q2, int *p)
{
    int ed[4], i, j, k, ts, dummy, el[2];
    el[0] = E;
    el[1] = F;
    k = 0;
    for (j = 0; j < 2; j++) {
        ed[0] = quad.elem[el[j]].frkt;
        ed[1] = quad.elem[el[j]].sckt;
        ed[2] = quad.elem[el[j]].trkt;
        ed[3] = quad.elem[el[j]].ftkt;
        for (i = 0; i < 4; i++)
            if ((ed[i] != q1) && (ed[i] != q2)) {
                ts = gonl_arra_govj(p, k, ed[i], &dummy);
                if (ts == 0) {
                    p[k] = ed[i];
                    k++;
                }
            }
    }
    if (k != 4) {
        fprintf(tmpout, "There must be four sur edges  [q1,q2]=[%d,%d]\n", q1, q2);
        for (i = 0; i < k; i++)
            fprintf(tmpout, "ed[%d]=%d\n", i, ed[i]);
        jans_disp_nudj(quad, E);
        jans_disp_nudj(quad, F);
        exit(0);
    }
}


void fupj_find_numk_jobd(efajor q, efajor * p)
{
    p->frvrt = q.frvrt;
    p->scvrt = q.scvrt;
    p->thvrt = q.thvrt;
    p->ftvrt = q.ftvrt;
    p->frkt = q.frkt;
    p->sckt = q.sckt;
    p->trkt = q.trkt;
    p->ftkt = q.ftkt;
}


void tolm_merg_wovm(fajor_sion3D * quad, int *corresp, int E, int F, float_curve * fc, int max_nb_smooth)
{
    int q1, q2, *p, a, b, c, d, i, j, k, obs, sd[2], rd[2];
    int loc[4], ts1, ts2, dummy, ned, *map_ed, e1, e2, e3, e4;
    int *map_qd, *map_nd, n1, n2, n3, n4, z;
    float_curve *temp_fc;
    kt *temp_ed;
    efajor *temp_qd;
    point *temp_nd;

    vunm_comm_dilp(*quad, E, F, &q1, &q2);
    p = (int *) malloc(4 * sizeof(int));
    dezq_surr_qong(*quad, E, F, q1, q2, p);
    sd[0] = quad->kt[q1].frvrt;
    sd[1] = quad->kt[q1].scvrt;
    rd[0] = quad->kt[q2].frvrt;
    rd[1] = quad->kt[q2].scvrt;
    obs = -1;
    for (i = 0; i < 2; i++) {
        for (j = 0; j < 2; j++)
            if (sd[i] == rd[j]) {
                obs = sd[i];
                break;
            }
        if (obs != -1)
            break;
    }
    if (obs == -1) {
        fprintf(tmpout, "Unable to find obs\n");
        exit(0);
    }
    loc[0] = quad->elem[E].frvrt;
    loc[1] = quad->elem[E].scvrt;
    loc[2] = quad->elem[E].thvrt;
    loc[3] = quad->elem[E].ftvrt;
    a = -1;
    b = -1;
    c = -1;
    d = -1;
    for (i = 0; i < 4; i++) {
        ts1 = gonl_arra_govj(sd, 2, loc[i], &dummy);
        ts2 = gonl_arra_govj(rd, 2, loc[i], &dummy);
        if ((ts1 == 0) && (ts2 == 0)) {
            a = loc[i];
            break;
        }
    }
    if (a == -1) {
        fprintf(tmpout, "Unable to find a\n");
        exit(0);
    }
    b = cuzh_next_kadm(*quad, E, a);
    c = cuzh_next_kadm(*quad, F, b);
    d = cuzh_next_kadm(*quad, F, c);

    quad->elem[E].frvrt = a;
    quad->elem[E].scvrt = b;
    quad->elem[E].thvrt = c;
    quad->elem[E].ftvrt = d;
    quad->elem[E].frkt = p[0];
    quad->elem[E].sckt = p[1];
    quad->elem[E].trkt = p[2];
    quad->elem[E].ftkt = p[3];
    for (i = 0; i < 4; i++) {
        if (quad->kt[p[i]].frent == F)
            quad->kt[p[i]].frent = E;
        if (quad->kt[p[i]].scent == F)
            quad->kt[p[i]].scent = E;
    }
    free(p);

    ned = quad->k_grs;
    temp_fc = (float_curve *) malloc(ned * sizeof(float_curve));
    for (i = 0; i < ned; i++)
        vewk_allo_jovk(max_nb_smooth, &temp_fc[i]);
    temp_ed = (kt *) malloc(ned * sizeof(kt));
    map_ed = (int *) malloc(ned * sizeof(int));
    k = 0;
    for (i = 0; i < ned; i++)
        if ((i != q1) && (i != q2)) {
            if (fc[i].st_grs >= max_nb_smooth) {
                fprintf(tmpout, "max_nb_smooth is exceeded\n");
                exit(0);
            }
            rekf_find_kujn_rukg(fc[i], &temp_fc[k]);
            temp_ed[k].frvrt = quad->kt[i].frvrt;
            temp_ed[k].scvrt = quad->kt[i].scvrt;
            temp_ed[k].frent = quad->kt[i].frent;
            temp_ed[k].scent = quad->kt[i].scent;
            map_ed[i] = k;
            k++;
        }
    map_ed[q1] = -1;
    map_ed[q2] = -1;

    for (i = 0; i < k; i++) {
        rekf_find_kujn_rukg(temp_fc[i], &fc[i]);
        quad->kt[i].frvrt = temp_ed[i].frvrt;
        quad->kt[i].scvrt = temp_ed[i].scvrt;
        quad->kt[i].frent = temp_ed[i].frent;
        quad->kt[i].scent = temp_ed[i].scent;
    }
    quad->k_grs = k;
    for (i = 0; i < ned; i++)
        lohm_dest_nosr(&temp_fc[i]);
    free(temp_fc);
    free(temp_ed);
    for (i = 0; i < quad->e_grs; i++) {
        e1 = quad->elem[i].frkt;
        quad->elem[i].frkt = map_ed[e1];
        e2 = quad->elem[i].sckt;
        quad->elem[i].sckt = map_ed[e2];
        e3 = quad->elem[i].trkt;
        quad->elem[i].trkt = map_ed[e3];
        e4 = quad->elem[i].ftkt;
        quad->elem[i].ftkt = map_ed[e4];
    }
    free(map_ed);

    temp_qd = (efajor *) malloc(quad->e_grs * sizeof(efajor));
    map_qd = (int *) malloc(quad->e_grs * sizeof(int));
    k = 0;
    for (i = 0; i < quad->e_grs; i++)
        if (i != F) {
            fupj_find_numk_jobd(quad->elem[i], &temp_qd[k]);
            map_qd[i] = k;
            k++;
        }

    for (i = 0; i < k; i++)
        fupj_find_numk_jobd(temp_qd[i], &quad->elem[i]);
    quad->e_grs = k;
    free(temp_qd);
    for (i = 0; i < quad->k_grs; i++) {
        e1 = quad->kt[i].frent;
        quad->kt[i].frent = map_qd[e1];
        e2 = quad->kt[i].scent;
        if (e2 != -1)
            quad->kt[i].scent = map_qd[e2];
    }
    free(map_qd);

    temp_nd = (point *) malloc(quad->n_grs * sizeof(point));
    map_nd = (int *) malloc(quad->n_grs * sizeof(int));
    k = 0;
    for (i = 0; i < quad->n_grs; i++)
        if (i != obs) {
            getf_find_rogc_todj(quad->knot[i], &temp_nd[k]);
            map_nd[i] = k;
            z = corresp[i];
            corresp[k] = z;
            k++;
        }

    for (i = 0; i < k; i++)
        getf_find_rogc_todj(temp_nd[i], &quad->knot[i]);
    quad->n_grs = k;
    free(temp_nd);
    for (i = 0; i < quad->e_grs; i++) {
        n1 = quad->elem[i].frvrt;
        quad->elem[i].frvrt = map_nd[n1];
        n2 = quad->elem[i].scvrt;
        quad->elem[i].scvrt = map_nd[n2];
        n3 = quad->elem[i].thvrt;
        quad->elem[i].thvrt = map_nd[n3];
        n4 = quad->elem[i].ftvrt;
        quad->elem[i].ftvrt = map_nd[n4];
    }
    for (i = 0; i < quad->k_grs; i++) {
        n1 = quad->kt[i].frvrt;
        quad->kt[i].frvrt = map_nd[n1];
        n2 = quad->kt[i].scvrt;
        quad->kt[i].scvrt = map_nd[n2];
    }
    free(map_nd);
}


void lujs_find_capm_nerv(fajor_sion3D * quad, int *corresp, float_curve * fc, int *nb_smooth, int max_nb_smooth)
{
    int e, E, F, ts, dummy, suc, i, nnd, nb_sm;
    nnd = quad->n_grs;
    nb_sm = *nb_smooth;
    fprintf(tmpout, "FIX MERGE\n");
    for (i = 0; i < nnd; i++) {
        suc = FAILURE;
        for (e = 0; e < quad->k_grs; e++)
            if (quad->kt[e].scent != -1) {
                E = quad->kt[e].frent;
                F = quad->kt[e].scent;
                ts = cojs_shar_sejm(*quad, E, F, &dummy);
                if (ts == 1) {
                    tolm_merg_wovm(quad, corresp, E, F, fc, max_nb_smooth);
                    suc = SUCCESS;
                    nb_sm = nb_sm - 2;
                }
            }
        if (suc == FAILURE)
            break;
    }
    *nb_smooth = nb_sm;
}
