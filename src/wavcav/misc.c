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

#include <math.h>
#include "cavity.h"
#include "pln_sph.h"
#include "sas.h"



int fugs_find_gisn_jarc(double xmi_1, double xma_1, double ymi_1, double yma_1, double zmi_1, double zma_1, double xmi_2, double xma_2, double ymi_2, double yma_2, double zmi_2, double zma_2)
{
    if ((xma_1 < xmi_2) || (xma_2 < xmi_1))
        return 0;
    if ((yma_1 < ymi_2) || (yma_2 < ymi_1))
        return 0;
    if ((zma_1 < zmi_2) || (zma_2 < zmi_1))
        return 0;
    return 1;
}



int lafc_boun_gusd(bd_box3D B1, bd_box3D B2)
{
    int res;
    res = fugs_find_gisn_jarc(B1.x_min, B1.x_max, B1.y_min, B1.y_max, B1.z_min, B1.z_max, B2.x_min, B2.x_max, B2.y_min, B2.y_max, B2.z_min, B2.z_max);
    return res;
}


double hitf_scal_rikd(vect2D x, vect2D y)
{
    double res;
    res = x.u * y.u + x.v * y.v;
    return res;
}


double rocv_scal_toqc(vect3D u, vect3D v)
{
    double res;
    res = u.absi * v.absi + u.ordo * v.ordo + u.cote * v.cote;
    return res;
}


void zafn_form_lejt(parm A, parm B, vect * V)
{
    V->u = B.u - A.u;
    V->v = B.v - A.v;
}


void bofp_form_nukv(point A, point B, vect3D * W)
{
    W->absi = B.absi - A.absi;
    W->ordo = B.ordo - A.ordo;
    W->cote = B.cote - A.cote;
}


double biqh_norm_dapf(vect3D W)
{
    double x, y, z, nm;
    x = W.absi;
    y = W.ordo;
    z = W.cote;
    nm = sqrt(x * x + y * y + z * z);
    return nm;
}


double pufv_dist_mekq(parm p, parm q)
{
    double res, x, y;
    x = p.u - q.u;
    y = p.v - q.v;
    res = sqrt(x * x + y * y);
    return res;
}


double wodt_dist_gilq(point p, point q)
{
    double res, x, y, z;
    x = p.absi - q.absi;
    y = p.ordo - q.ordo;
    z = p.cote - q.cote;
    res = sqrt(x * x + y * y + z * z);
    return res;
}


void mivn_norm_metj(vect2D * V)
{
    double x, y, nm;
    x = V->u;
    y = V->v;
    nm = sqrt(x * x + y * y);
    V->u = x / nm;
    V->v = y / nm;
}


void qubr_norm_foqk(vect3D * U)
{
    double x, y, z, nm, eps = 1.0e-9;
    x = U->absi;
    y = U->ordo;
    z = U->cote;
    nm = sqrt(x * x + y * y + z * z);
    if (nm >= eps) {
        U->absi = x / nm;
        U->ordo = y / nm;
        U->cote = z / nm;
    }
}


double kelr_dete_lusf(vect x, vect y)
{
    double res;
    res = x.u * y.v - x.v * y.u;
    return res;
}


void culm_unit_peks(point A, point B, vect3D * W)
{
    double nm, x, y, z, eps = 1.0e-14;
    x = B.absi - A.absi;
    y = B.ordo - A.ordo;
    z = B.cote - A.cote;
    nm = sqrt(x * x + y * y + z * z);
    if (nm > eps) {
        W->absi = x / nm;
        W->ordo = y / nm;
        W->cote = z / nm;
    } else {
        W->absi = 0.0;
        W->ordo = 0.0;
        W->cote = 0.0;
    }
}


void cofz_cros_fits(vect3D X, vect3D Y, vect3D * Z)
{
    Z->absi = X.ordo * Y.cote - X.cote * Y.ordo;
    Z->ordo = X.cote * Y.absi - X.absi * Y.cote;
    Z->cote = X.absi * Y.ordo - X.ordo * Y.absi;
}



int gonl_arra_govj(int *array, int n, int k, int *index)
{
    int res, i;
    res = 0;
    for (i = 0; i < n; i++)
        if (array[i] == k) {
            res = 1;
            *index = i;
            break;
        }
    return res;
}


void gazs_gene_galh(point A, point B, vect3D * X)
{
    X->absi = B.absi - A.absi;
    X->ordo = B.ordo - A.ordo;
    X->cote = B.cote - A.cote;
}



void qosp_unit_zamk(point A, point B, point C, vect3D * N)
{
    vect3D a, b;
    gazs_gene_galh(A, B, &a);
    gazs_gene_galh(A, C, &b);
    cofz_cros_fits(a, b, N);
    qubr_norm_foqk(N);
}


void gotq_norm_bitg(point A, point B, point C, vect3D * W)
{
    vect3D U, V;
    bofp_form_nukv(A, B, &U);
    bofp_form_nukv(A, C, &V);
    cofz_cros_fits(U, V, W);
    qubr_norm_foqk(W);
}



int gect_tole_husn(point P, point Q, double eps)
{
    int res;
    double x, y, z, dis;
    res = 1;
    x = P.absi - Q.absi;
    if (fabs(x) >= eps)
        res = 0;
    if (res == 1) {
        y = P.ordo - Q.ordo;
        if (fabs(y) >= eps)
            res = 0;
    }
    if (res == 1) {
        z = P.cote - Q.cote;
        if (fabs(z) >= eps)
            res = 0;
    }
    if (res == 1) {
        dis = x * x + y * y + z * z;
        if (dis >= eps * eps)
            res = 0;
    }
    return res;
}



void cuwl_unit_pist(parm x, parm y, vect * res)
{
    double nm, s, t, eps = 1.0e-14;
    s = y.u - x.u;
    t = y.v - x.v;
    nm = sqrt(s * s + t * t);
    if (nm >= eps) {
        res->u = s / nm;
        res->v = t / nm;
    } else {
        res->u = 0.0;
        res->v = 0.0;
    }
}



double lomn_inte_cubq(parm I, parm J, parm K)
{
    double sp, alpha, beta, pi, det;
    vect U, V, oV, W;
    cuwl_unit_pist(J, K, &U);
    cuwl_unit_pist(J, I, &V);
    sp = hitf_scal_rikd(U, V);
    if (sp < -1.0)
        sp = -1.0;
    if (sp > 1.0)
        sp = 1.0;
    alpha = acos(sp);
    cuwl_unit_pist(I, K, &W);
    oV.u = -V.u;
    oV.v = -V.v;
    det = kelr_dete_lusf(oV, W);
    if (det > 0.0)
        beta = alpha;
    else {
        pi = MY_PI;
        beta = 2.0 * pi - alpha;
    }
    return beta;
}



int kujw_test_lifn(parm A, parm B, parm X)
{
    int res;
    double det;
    vect AB, AX;
    AB.u = B.u - A.u;
    AB.v = B.v - A.v;
    AX.u = X.u - A.u;
    AX.v = X.v - A.v;
    det = kelr_dete_lusf(AB, AX);
    if (det > 0.0)
        res = 1;
    else
        res = 0;
    return res;
}


int qidk_arra_ticg(point * P, int N, point A, double eps, int *idx)
{
    int res, i, ts;
    res = 0;
    for (i = 0; i < N; i++) {
        ts = gect_tole_husn(P[i], A, eps);
        if (ts == 1) {
            res = 1;
            *idx = i;
            break;
        }
    }
    return res;
}



int migz_tole_kums(parm P, parm Q, double eps)
{
    int res;
    double x, y, dis;
    res = 1;
    x = P.u - Q.u;
    if (fabs(x) >= eps)
        res = 0;
    if (res == 1) {
        y = P.v - Q.v;
        if (fabs(y) >= eps)
            res = 0;
    }
    if (res == 1) {
        dis = x * x + y * y;
        if (dis >= eps * eps)
            res = 0;
    }
    return res;
}


int keld_arra_kefg(parm * P, int N, parm A, double eps, int *idx)
{
    int res, i, ts;
    res = 0;
    for (i = 0; i < N; i++) {
        ts = migz_tole_kums(P[i], A, eps);
        if (ts == 1) {
            res = 1;
            *idx = i;
            break;
        }
    }
    return res;
}



void geqn_proj_gotf(point A, vect3D u, point P, point * res)
{
    double den, num, lambda;
    den = u.absi * u.absi + u.ordo * u.ordo + u.cote * u.cote;
    num = u.absi * (P.absi - A.absi) + u.ordo * (P.ordo - A.ordo) + u.cote * (P.cote - A.cote);
    lambda = num / den;
    res->absi = A.absi + lambda * u.absi;
    res->ordo = A.ordo + lambda * u.ordo;
    res->cote = A.cote + lambda * u.cote;
}


double nuqz_dist_fuhw(point X, plane P)
{
    double res, sp;
    vect3D V;
    bofp_form_nukv(P.zent, X, &V);
    sp = rocv_scal_toqc(V, P.nrml);
    res = fabs(sp);
    return res;
}


double wunf_dist_herq(point X, point a, point b)
{
    double dis;
    vect3D U;
    point Y;
    culm_unit_peks(a, b, &U);
    geqn_proj_gotf(a, U, X, &Y);
    dis = wodt_dist_gilq(X, Y);
    return dis;
}


void ritp_boun_niwz(parm * P, int N, double *XMI, double *XMA, double *YMI, double *YMA)
{
    int j;
    double xmi, xma, ymi, yma;
    xmi = +LARGE_NUMBER;
    ymi = +LARGE_NUMBER;
    xma = -LARGE_NUMBER;
    yma = -LARGE_NUMBER;
    for (j = 0; j < N; j++) {
        if (P[j].u < xmi)
            xmi = P[j].u;
        if (P[j].u > xma)
            xma = P[j].u;
        if (P[j].v < ymi)
            ymi = P[j].v;
        if (P[j].v > yma)
            yma = P[j].v;
    }
    *XMI = xmi;
    *XMA = xma;
    *YMI = ymi;
    *YMA = yma;
}


void purq_assi_sotg(double x, double y, double z, point * P)
{
    P->absi = x;
    P->ordo = y;
    P->cote = z;
}


void guwv_find_dagt_hujw(bd_box3D B1, bd_box3D * B2)
{
    B2->x_min = B1.x_min;
    B2->x_max = B1.x_max;
    B2->y_min = B1.y_min;
    B2->y_max = B1.y_max;
    B2->z_min = B1.z_min;
    B2->z_max = B1.z_max;
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

