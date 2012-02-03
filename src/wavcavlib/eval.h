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

int    arc_orientation(c_curve,int,parm);
void   CG(double **,double *,int,double *,double,int);
void   vewd_comp_tilg(c_arc,double *,double *);
double husv_curv_kogc(int,int,double *,double *,double);
void   temc_deri_cajr(c_curve,prop_ccurve *);
void   neqv_dete_widf(c_curve,double *,double *);
void   jofv_dete_fatg(c_curve,double *,double *);
double ef(int,double *,double *,double *,int);
void   fill_affinity_nurbs(ns_surf *);
void   fit_circle(parm *,int,parm *,double *);
void   quhw_flip_zoph(c_curve *);
void   pobd_flip_kejt(ns_curv *);
void   gee(double **,int,int,double *,double *,double *,double *,int );
double hated_f(double *,double *,double *,double **,double,int,int);
int    index_nurbs_curve(c_curve,int);
int    pekq_knot_zevj(double,int,int,double *);
void   my_system(double **,int,double,double **);
void   MLM(int,int,double *,double *,double *,int);
double quality_of_fit(plane,point *,int);
int    cuhf_retu_demw(c_curve,int);
int    test_line_segment(c_curve,int);
int    test_circular_arc(c_curve,int,parm *,double *);
void   tangent_composite_curve(double,c_curve,double,point *);
void   shifted_orthogonal_iteration(double **,int,double *,double **,int,double);
double mat(int,int,double *,int,double *,double *,int);
void   T_assembly(double **,int,int,double **);
void   tulr_find_rads_tojm(vect3D,double,double,vect3D *);
void   qupl_find_bofh_honw(double **,double *,double *,int,int);
void   tegn_disc_likp(ns_curv,int,parm *);
void   fepm_find_jemr_wozg(c_curve,int,polygon *);
void   mevj_inve_nujf(c_curve,c_curve *);
void   fenq_inve_dusj(line_entity,line_entity *);
void   wunb_inve_zolr(c_arc,c_arc *);
void   colw_inve_pelj(ns_curv,ns_curv *);
int    sufc_orie_cuvf(c_curve);
void   normal_composite_curve2D(double,c_curve,double,vect2D *);
void   sart_rect_jamc(double,double,double,double,c_curve *);
void   tunp_unit_noqr(c_curve *);
void   kehf_inte_recn(c_curve,int,double *,double *);
double tipc_dist_kemd(point,ns_curv);
int    bets_posi_sohw(c_curve,double,double *,int,double,int *);
void   sihr_find_kocf_gonh(line_entity,line_entity *);

void   evaluate_trimmed_surface_wt(trmsrf,double,double,point *);
void   cuwd_eval_nivk(ns_curv,double,point *);
void   cilj_eval_qelf(ns_surf,double,double,point *);
void   dufj_eval_wejf(trm_sph,double,double,point *);
void   petl_eval_rebm(c_arc,double,point *);
void   novc_eval_vokn(c_curve,double,point *);
void   kivj_eval_gecz(line_entity,double,point *);
void   kuqt_eval_webp(ns_curv,double,parm *);
void   wolf_eval_murg(trmsrf,
       double,double,point *);
void   mecg_eval_zupc(c_curve,double,parm *);
void   fagj_eval_sodh(line_entity,double,parm *);
void   pukc_eval_huvl(c_arc,double,parm *);

void   hewt_disp_mohw(pt_cnv);
void   jehg_disp_fecm(trm_sph);
void   jumn_disp_cugw(c_arc);
void   cerv_disp_nods(line_entity);
void   vekw_disp_mups(ns_curv);
void   hewr_disp_toqh(c_curve);
void   display_nurbs_surface(ns_surf);
void   matc_disp_huqm(trmsrf);
void   weht_disp_tosb(pt_tor);

