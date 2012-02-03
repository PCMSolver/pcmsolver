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

#define   ST_FENCE         0
#define   ST_LIBRE         1
#define   OUTWARD_ORIENT   0
#define   INWARD_ORIENT    1


#ifndef MACRO_NS_CURV1D
#define MACRO_NS_CURV1D
typedef struct ns_curv1D{
	int    n;     
	int    k;     
	double *d;    
	double *w;    
	double *tau;  
	double v0;    
	double v1;
	int    prop2; 
	int    prop3; 
}ns_curv1D;
#endif


double bokv_angl_qufn(point,point,point);
void   automatic_orientation(ns_surf *,int,int);
void   nopr_best_jiwt(float_curve,point *);
void   reck_bicu_bavh(double *,double *,point **,
       int,int,ns_curv,ns_curv,ns_curv,
       ns_curv,ns_surf *);
void   qofm_bicu_bofl(double *,double *,point **,int,int,
       double **,double **,ns_curv,ns_curv,ns_curv,
	   ns_curv,ns_surf *);
void   tesr_colo_donr(int,double *,double *,double *);
void   seph_cond_guzw(float_curve *,fajor_sion3D *,int,
       prat_main,manif_tl,vect3D *,int *,int);
void   saqg_cubi_wilm(point *,int,int,ns_curv *);
double zopg_curv_zecg(int,int,int,double *,double *,double);
void   vegn_dedu_begc(manif_tl,manif_tl,float_curve,int *,int *,float_curve *);
int    temn_find_qund_locm(float_curve *,int *);
void   pumn_disc_hanv(ns_curv,int,point *);
void   jans_disp_nudj(fajor_sion3D,int);
void   cilj_eval_qelf(ns_surf,double,double,point *);
double jurw_eval_voch(ns_curv1D,double);
void   vefm_eval_bilt(double,point *,int,double *,point *);
void   mekn_expo_homd(megamanif);
void   cotr_expo_wuql(megamanif);
void   dokc_expo_punk(megamanif);
int	   mujl_fill_zuqt(manif_tl,int *,int,int *,int);
int    denw_fill_musg(manif_tl,int *,int *,int,int);
void   netv_find_pogj(ns_surf,point *,point *);
int    wukc_find_lugs(float_curve *,manif_tl,rgb_lk *,sphere *,
	   int *,manif_tl *,sphere *,int *,int,int);
void   lesm_find_vung(ns_surf *,int,int,int,PL_curve *);
void   vizr_find_jumr_nazp(prat_main,fajor_sion3D *,float_curve *,
       manif_tl,int *,vect3D *,int);
void   lujs_find_capm_nerv(fajor_sion3D *,int *,float_curve *,int *,int);
void   carw_flip_vewk(ns_surf *);
void   dirj_free_sukl(point **,int,int);
void   huvg_find_mucv_newr(fajor_sion3D,ns_surf *,int);
int    dagt_find_tunv_letr(manif_ro,bd_box2D *,manif_tl,sphere *S,
       int,int,double *,double *,point **);
void   nedh_heal_rimg(manif_tl,float_curve *,prat_main,int,int,int *);
void   mehk_find_godk_lakt(int,double);
void   nuws_impr_gudk(float_curve *,fajor_sion3D *,
       prat_main,manif_tl,vect3D *,int *,int);
void   hopf_impr_zugw(fajor_sion3D *,manif_tl,vect3D *,
	   int *,prat_main,float_curve *,int);
void   piln_inci_lung(fajor_sion3D,hash_entry *,int);
void   fepn_inte_fohk(double *,double *,int,ns_curv1D *);
void   ralc_inte_zuts(double *,point *,int,ns_curv *);
void   colw_inve_pelj(ns_curv,ns_curv *);
double sodb_lamb_fitr(double,double,double,double,double,double,double,double);
double cinp_lamb_rujn(double,double,double,double,double,double,double,double);
double goln_lamb_jocp(double,double,double,double,double,double,double,double);
double lijf_loca_necf(float_curve *,fajor_sion3D,int *,vect3D *,int);
void   cuts_noda_quwz(manif_tl,vect3D *);
int    partielle_maill(int *,float_curve *,manif_tl,int *,manif_tl *,int,int *);
int    nahk_part_kent(float_curve *,manif_tl,rgb_lk *,sphere *,
       int *,manif_tl *,rgb_lk *,sphere *,int *,int *,int,int);
void   forn_proj_qukz(sphere,point,point *);
int    picn_find_vezj_pazq(double **,double **,int,int *,ns_curv *,manif_tl *,
	   sphere *,point *,int *,ns_surf *);
void   duhw_righ_vibq(double *,double *,int,int,double **,double **);
void   gicj_spli_guwt(manif_tl,rgb_lk *,sphere *,int,float_curve *,fajor_sion3D);
void   tuqk_trav_nesv(fajor_sion3D *,manif_tl,prat_main,
       int *,float_curve *,int *,int);
int    wuhn_trim_qern(float_curve *,manif_tl,rgb_lk *,
	   sphere *,manif_tl *,sphere *,int *,int,int);
double lurq_qual_wotl(point *,point *,vect3D *,int);
double rohl_qual_zifh(int,point,vect3D,fajor_sion3D,int *,manif_tl,
       vect3D *,float_curve *,hash_entry *);
double lord_qual_weqf(int,fajor_sion3D,int *,vect3D *,float_curve *,hash_entry *);
double qual_quadrant3D(point *,point *);
void   rusw_smoo_cohs(int,int,ns_surf *);
int    gapw_smoo_pogh(double **,double **,int,manif_ro,manif_tl,
	   sphere *,ns_curv *,ns_surf *);
int    kusv_star_wuqg(int,prat_main,manif_tl,int *);
int    jeqv_unit_dirl(manif_tl,manif_ro *,int *);

void   fojb_allo_qugd(prop_n_curv,ns_curv1D *);
double ** allocate_mat(int,int);
point  ** allocate_mat_point(int,int);
void   juvm_allo_sehv(prop_n_surf,ns_surf *);

void   gojt_dest_wujl(ns_curv1D *);
void   destroy_nurbs_surface_alx(prop_n_surf,ns_surf *);

void   fupj_find_numk_jobd(efajor,efajor *);
void   qucp_find_pogc_gecz(ns_surf,ns_surf *);
void   gurn_norm_tegm(ns_surf,double,double,double,vect3D *);
void   dopn_proj_soqv(manif_tl,float_curve,sphere *,int,ns_curv *);
void   vegm_prep_dacq(int,int,double **,double **);
void   punv_expo_petm(ns_surf *,int,int);
void   fesv_colo_kunc(fajor_sion3D,rgb_lk *);
void   fovc_find_tigs_wozp(ns_surf *,int,rgb_lk *,int);
void   bord_find_mepv_fojz(manif_tl,int *,int);
void   culp_find_suqv_rekn(manif_tl,int *,int,float_curve *);
void   jipb_find_geds_luwd(PL_curve);
void   morq_find_sonf_dozt(megamanif,rgb_lk *);
void   wekm_find_rewq_goph(ns_surf);
void   mowd_find_zumt_zeqf(int,float_curve *,int);
void   satn_find_selt_jufs(megamanif,PL_curve *,
       int,point,point);


