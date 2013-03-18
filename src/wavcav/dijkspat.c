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
#include "cavity.h"
#include "splinemol.h"
#include "geodesic.h"


int monh_extr_jaws(int *Q, int lQ, double *d)
{
    int u, i, r, s;
    double sm = 1.0e+10;

    for (i = 0; i < lQ; i++) {
        r = Q[i];
        if (d[r] < sm) {
            sm = d[r];
            u = r;
        }
    }

    for (i = 0; i < lQ; i++)
        if (Q[i] == u) {
            s = i;
            break;
        }
    for (i = s; i < lQ; i++)
        Q[i] = Q[i + 1];
    return u;
}


int kegj_test_cigk(int *x, int l, int r)
{
    int ts = 0, i;
    for (i = 0; i < l; i++)
        if (x[i] == r) {
            ts = 1;
            break;
        }
    return ts;
}


int fips_rela_gusl(int u, int *pi, double *d, int *S, int lS, int *Q, int lQ, prat_main g)
{
    int i, j, n, e, v, ts, res, c, pr;
    n = g.dgr[u];
    c = lQ;
    res = 0;
    for (i = 0; i < n; i++) {
        e = g.incd[u][i];
        if (g.kt[e].str == u)
            v = g.kt[e].ter;
        else
            v = g.kt[e].str;

        ts = 0;
        for (j = 0; j < lS; j++)
            if (S[j] == v) {
                ts = 1;
                break;
            }
        if (ts == 0) {
            if (d[v] > d[u] + g.gew[e]) {
                d[v] = d[u] + g.gew[e];
                pi[v] = u;
                pr = kegj_test_cigk(Q, c, v);
                if (pr == 0) {
                    Q[c] = v;
                    c++;
                    res++;
                }
            }
        }
    }
    return res;
}



void cofd_dijk_nugc(prat_main g, int s, double *d, int *pi)
{
    int i, n, *S, *Q, lS, lQ, u, ts, emptyness;
    double infinity = 1.0e+4;

    n = g.v_grs;
    for (i = 0; i < n; i++)
        d[i] = infinity;
    lS = 0;
    lQ = 0;
    Q = (int *) malloc(2 * n * sizeof(int));
    S = (int *) malloc(2 * n * sizeof(int));
    Q[0] = s;
    lQ = 1;
    d[s] = 0.0;
    emptyness = 0;
    while (emptyness == 0) {
        u = monh_extr_jaws(Q, lQ, d);
        lQ--;
        S[lS] = u;
        lS++;
        ts = fips_rela_gusl(u, pi, d, S, lS, Q, lQ, g);
        lQ = lQ + ts;
        if (lQ == 0)
            emptyness = 1;
    }
    free(Q);
    free(S);
    pi[s] = s;
}



int qofj_find_wogv(prat_main g, int nd1, int nd2)
{
    int e = -1, i, n;
    n = g.dgr[nd1];
    for (i = 0; i < n; i++) {
        e = g.incd[nd1][i];
        if ((g.kt[e].str == nd2) || (g.kt[e].ter == nd2))
            break;
    }
    return e;
}



void jalf_gene_homz(prat_main g, int source, int dest, int *station, int *l)
{
    int *pi, n, i, k, *temp, L;
    double *d;
    n = g.v_grs;
    pi = (int *) malloc(n * sizeof(int));
    d = (double *) malloc(n * sizeof(double));
    cofd_dijk_nugc(g, source, d, pi);
    temp = (int *) malloc(n * sizeof(int));
    temp[0] = dest;
    for (i = 0; i < n; i++) {
        if (temp[i] == source) {
            L = i + 1;
            break;
        }
        k = pi[temp[i]];
        temp[i + 1] = k;
    }
    free(pi);
    free(d);
    for (i = 0; i < L; i++)
        station[i] = temp[L - i - 1];
    free(temp);
    *l = L;
}
