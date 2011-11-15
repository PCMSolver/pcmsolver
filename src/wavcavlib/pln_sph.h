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

#define  INFINITE_VOR_VER    0    
#define  MAX_MEMBERS         100  
#define  MAX_ARCS            75   
#define  MAX_TRIM_SURF1      1750 
#define  MAX_TRIM_SURF2      7000 
#define  MAX_TRIM_SURF3      8500 
#define  MAX_TET_INC         100  
#define  MAXCOMP             100   
#define  MAX_INTERNAL_CURVES 6    
#define  MAX_ATOMS           3800  
#define  IGS_MAX_NND         5000 
#define  IGS_MAX_NEL         3500 
#define  IGS_MAX_NED         4500 
#define  MAX_BLEND_TOR       6500 
#define  N_COMPLETE_CIR      8    
#define  MAX_LOC_BLEND       175  
#define  MAX_SELFINTERS      1000
#define  MAX_TOR_PATCH       4000
#define  VERBOSE             1
#define  SILENT				 0
#define  CLOCKWISE           1
#define  COUNTER_CLOCKWISE   2
#define  SUCCESS             1
#define  FAILURE             0
#define  STRICTLY_INSIDE     0
#define  STRICTLY_OUTSIDE    1
#define  ON_BOUNDARY         2
#define  INTERNAL_VOR_VER    1  
 

int verbose_variable;

 
#ifndef MACRO_PLANE
#define MACRO_PLANE
typedef struct plane{
	vect3D  nrml;	
	point   zent; 
}plane;
#endif



#ifndef MACRO_C_ARC3D
#define MACRO_C_ARC3D
typedef struct c_arc3D{
	int    c_cir;
	point  zent;      
	double rad;      
	vect3D nrml;      
	point  begn;       
	point  term;        
}c_arc3D; 
#endif


#ifndef MACRO_CIRCLE2D
#define MACRO_CIRCLE2D
typedef struct circle2D{
	parm   zent;
	double rad;
}circle2D; 
#endif


#ifndef MACRO_KT_TRI
#define MACRO_KT_TRI
typedef struct kt_tri{
	int st_node;   
	int tr_node;   
}kt_tri,arrow;   
#endif


#ifndef MACRO_PL_CURVE2D
#define MACRO_PL_CURVE2D
typedef struct PL_curve2D{
	int  v_grs;
	parm *vertex;
}PL_curve2D;
#endif

 
#ifndef MACRO_SET_ARCS
#define MACRO_SET_ARCS
typedef struct set_arcs{
	int            ar_grs;
	c_arc3D *C;
	int            *par_idx;     
}set_arcs;
#endif


#ifndef MACRO_PT_TOR
#define MACRO_PT_TOR
typedef struct pt_tor{
	c_arc3D alpha;
	c_arc3D beta;
	c_arc3D gamma;
	c_arc3D delta;
}pt_tor;
#endif


#ifndef MACRO_TRM_SPH
#define MACRO_TRM_SPH
typedef struct trm_sph{
	point  zent;
	double rad;
	double beta;  
	vect3D nrml;
}trm_sph; 
#endif


#ifndef MACRO_ADJ_ENTRY
#define MACRO_ADJ_ENTRY
typedef struct adj_entry{
	int nb_neighbors;
	int *neighbor;
}adj_entry;
#endif


#ifndef MACRO_ADJ_HASH
#define MACRO_ADJ_HASH
typedef struct adj_hash{
	int             nb_spheres;        
	adj_entry *entry;    
	adj_entry *inter;    
}adj_hash;
#endif


#ifndef MACRO_LINE_ENTITY
#define MACRO_LINE_ENTITY
typedef struct line_entity{
	point p1;
	point p2;
}line_entity;
#endif



#ifndef MACRO_XFORM
#define MACRO_XFORM
typedef struct xform{
	double R[3][3];
	double T[3];
}xform;
#endif



#ifndef MACRO_C_ARC
#define MACRO_C_ARC
typedef struct c_arc{
	double ZT;  
	point  p1;  
	point  p2;  
	point  p3;  
	double A;
	double B;
	xform  TR;  
}c_arc;
#endif



#ifndef MACRO_C_CURVE
#define MACRO_C_CURVE
typedef struct c_curve{
	int          N;     
	int          *type; 
	int          nle;   
	int          nca;   
	int          nnc;   
	line_entity  *le;   
	c_arc *ca;   
	ns_curv  *nc;   
}c_curve;
#endif


#ifndef MACRO_PT_CNV
#define MACRO_PT_CNV
typedef struct pt_cnv{
	int             sitn;    
	double          scl;     
	c_arc3D  alpha;
	c_arc3D  beta;
	c_arc3D  gamma;
	trm_sph  T;   
	c_curve cc; 
}pt_cnv;
#endif


#ifndef MACRO_JOINTBEZ
#define MACRO_JOINTBEZ
typedef struct jointbez{
	bez_crv B_beg;
	bez_crv B_end;
}jointbez;
#endif



#ifndef MACRO_BEZ_CRV2D
#define MACRO_BEZ_CRV2D
typedef struct bez_crv2D{
	int   dgr;
	parm  *ctr;
}bez_crv2D;
#endif


#ifndef MACRO_JOINTBEZ2D
#define MACRO_JOINTBEZ2D
typedef struct jointbez2D{
	bez_crv2D B_beg;
	bez_crv2D B_end;
}jointbez2D;
#endif


#ifndef MACRO_BLEND_NONSELF
#define MACRO_BLEND_NONSELF
typedef struct blend_nonself{
	jointbez	   alpha;
	c_arc3D beta;
	jointbez	   gamma;
	c_arc3D delta;
}blend_nonself;
#endif
  


#ifndef MACRO_TRMSRF
#define MACRO_TRMSRF
typedef struct trmsrf{
	int             type;      
	int             boundary;  
	int             nb_inner;  
	trm_sph  ts;        
	pt_tor  pt;        
	pt_cnv   pc;        
	blend_nonself   bn;        
	c_curve cc;        
	c_curve *inner;    
	c_curve old_cc;    
}trmsrf;
#endif


#ifndef MACRO_PROP_CCURVE
#define MACRO_PROP_CCURVE
typedef struct prop_ccurve{
	int N;          
	int nle;        
	int nca;        
	int nnc;        
	int n[MAXCOMP]; 
	int k[MAXCOMP]; 
}prop_ccurve;
#endif



#ifndef MACRO_TRM_IDS
#define MACRO_TRM_IDS
typedef struct trm_ids{
	int nb_inter;	       
	int *nb_curve_comp;    
	int **list_curve_comp; 
}trm_ids;
#endif


#ifndef MACRO_RAT_BEZ
#define MACRO_RAT_BEZ
typedef struct rat_bez{
	int    n;      
	point  *cont;  
	double *w;     
}rat_bez;
#endif


#ifndef MACRO_POLYGON
#define MACRO_POLYGON
typedef struct polygon{
	int   v_grs;          
	parm  *vertex;              
	int   nb_inner_boundaries;  
	int   *nb_local_vertices;   
}polygon;
#endif



#ifndef  MACRO_MULT_CONN
#define  MACRO_MULT_CONN
typedef struct mult_conn{
	int    v_grs;      
	parm   *vertex;		     
	double *zt;            
	int    *flag;            
	int    *mapglob;         
	int    nb_vr_outer;      
	int    nb_inner_polygons;
	int    *nb_vr_inner;     
}mult_conn;
#endif



#ifndef  MACRO_PAIR_NODES
#define  MACRO_PAIR_NODES
typedef struct pair_nodes{
	int nb_pairs;       
	int *first;         
	int *second;        
	
	int *type_cut;      
	int *supp_sphere;   
}pair_nodes;
#endif


#ifndef  MACRO_PLG3D
#define  MACRO_PLG3D
typedef struct plg3D{
	int v_grs;           
	point *vertex;              
	int nb_inner_boundaries;   
	int *nb_local_vertices;    
}plg3D;
#endif


#ifndef  MACRO_BDR3D
#define  MACRO_BDR3D
typedef struct bdr3D{
	plg3D interp;		
	plg3D offset_in;    
	plg3D offset_out;   
}bdr3D;
#endif


#ifndef MACRO_CLD_PLG
#define MACRO_CLD_PLG
typedef struct cld_plg{
	int  v_grs;	
	parm *vertex;       
	int  *type;			
}cld_plg;
#endif


#ifndef MACRO_BARYC2D
#define MACRO_BARYC2D
typedef struct baryc2D{
	double lambda1;
	double lambda2;
	double lambda3;
}baryc2D;
#endif


#ifndef MACRO_JOINTBEZ
#define MACRO_JOINTBEZ
typedef struct jointbez{
	bez_crv B_beg;
	bez_crv B_end;
}jointbez;
#endif

  

#ifndef MACRO_BEZ_CRV2D
#define MACRO_BEZ_CRV2D
typedef struct bez_crv2D{
	int   dgr;
	parm *ctr;
}bez_crv2D;
#endif


#ifndef MACRO_JOINTBEZ2D
#define MACRO_JOINTBEZ2D
typedef struct jointbez2D{
	bez_crv2D B_beg;
	bez_crv2D B_end;
}jointbez2D;
#endif


#ifndef MACRO_BLEND_NONSELF
#define MACRO_BLEND_NONSELF
typedef struct blend_nonself{
	jointbez	   alpha;
	c_arc3D beta;
	jointbez	   gamma;
	c_arc3D delta;
}blend_nonself;
#endif


#ifndef MACRO_FAST_ARC_COMP
#define MACRO_FAST_ARC_COMP
typedef struct fast_arc_comp{
	double phi,theta;
	double alpha_s,alpha_t;
	point mid;
	bd_box3D B;
}fast_arc_comp;
#endif


int    keld_arra_kefg(parm *,int,parm,double,int *);
int    gonl_arra_govj(int *,int,int,int *);
int    rijc_find_podl_zigj(c_arc3D,point);
int    ruqs_find_wazv_detp(c_arc3D,fast_arc_comp,point,double);
int    qidk_arra_ticg(point *,int,point,double,int *);
void   purq_assi_sotg(double,double,double,point *);
double begj_blen_nugz(double);
void   ritp_boun_niwz(parm *,int,double *,double *,double *,double *);
void   homs_boun_gosm(point *,int,double *,double *,
	   double *,double *,double *,double *);
void   bemh_find_nuwl_bafj(PL_curve *,int,double *,double *,
       double *,double *,double *,double *);
void   CG(double **,double *,int,double *,double,int);
int    riwf_chec_zibc(polygon,double,parm *);
int    sulw_circ_mojt(parm,double,parm,double,parm *);
void   cehk_circ_jesw(parm ,parm ,parm,double *,parm *);
void   tesr_colo_donr(int,double *,double *,double *);
void   compute_normal(point,point,point,vect3D *);
void   cofz_cros_fits(vect3D,vect3D,vect3D *);
double kelr_dete_lusf(vect,vect);
void   pork_find_nogk_qijr(c_arc3D,fast_arc_comp *);
double hefk_find_cejs_vuhl(int);
void   mopb_dire_woqp(vect3D,double,double,vect3D *);
void   medj_disc_mupq(circle3D,int,PL_curve *);
void   vuch_disc_mogv(c_arc3D,int,PL_curve *);
void   display_plane(plane);
void   modg_disp_huwt(manif_tl);
void   display_circle3D(circle3D);
void   nepf_disp_bulp(c_arc3D);
double cijv_dist_laph(point,c_arc3D,int *);
double cutj_dist_rulb(point,point,c_arc3D);
double fozm_dist_lojn(c_arc3D,c_arc3D);
double regv_dist_nomw(c_arc3D,c_arc3D);
double pufv_dist_mekq(parm,parm);
double wodt_dist_gilq(point,point);
double wunf_dist_herq(point,point,point);
double nuqz_dist_fuhw(point,plane);
void   wusd_eval_jomk(c_arc3D,double,double,double,point *);
void   cutn_eval_mecn(blend_nonself,parm,point *);
void   nehl_eval_segt(c_arc3D,double,point *);
void   rows_eval_qusg(pt_cnv,parm,point *);
void   rogv_eval_dukw(pt_tor,parm,point *);
void   vumf_find_pobs_gafp(manif_tl,int);
void   rinf_find_qijc_wuvr(manif_tl,int);
void   jegn_find_wikc_poqc(manif_tl,int);
void   export_cloud_voronoi(point *,int);
void   josh_expo_qazg(point,double);
int    extr_voronoi_vertices(point *);
void   zafn_form_lejt(parm,parm,vect *);
void   bofp_form_nukv(point,point,vect3D *);
void   tehg_free_dacp(double **,int,int);
int    javs_find_gokm_pedt(double **,double *,double *,int);
void   gazs_gene_galh(point,point,vect3D *);
void   hetc_glob_kaqc(double,double,sphere *,int,adj_hash,int *);
void   tojw_hash_pojh(double,sphere *,int,adj_hash *,int,int);
void   integrity_offset_dir(manif_tl,plane *,int);
void   integrity_offset_dir_once(manif_tl,plane);
double lomn_inte_cubq(parm,parm,parm);
int    nefr_inte_sujc(point *,double,sphere *,sphere *);
void   laqv_inte_jotl(c_arc2D,double *,double *);
void   qirp_inte_ligr(c_arc3D,double *,double *);
void   fesg_inve_pahj(vect3D,double,double,vect3D *);
double vupq_lamb_qofc(double,double,double,double,double,double,double,double);
double dopg_lamb_nupd(double,double,double,double,double,double,double,double);
double mofr_lamb_powg(double,double,double,double,double,double,double,double);
void   reds_load_katr(char *,int,sphere *);
double tepc_leng_ziql(c_arc3D);
double cuhn_find_rewp_cefm(double *,int);
void   mowb_mesh_wefs(point,double,int,int,int,manif_tl *);
void   wogn_mesh_tuhl(point,double,int,int,manif_tl *);
void   renw_midp_mocw(c_arc3D,point *);
void   demq_mole_lufh();
void   lujg_find_pikz_figc(int,int);
void   gotq_norm_bitg(point,point,point,vect3D *);
double biqh_norm_dapf(vect3D);
void   mivn_norm_metj(vect2D *);
void   qubr_norm_foqk(vect3D *);
int    vofr_numb_qimf(char *, atom **);
int    sogf_para_lekj(double,sphere,sphere,circle3D *,circle3D *);
void   qezj_find_tukd_wesg(prop_ccurve *);
int    bosr_plan_revg(point,vect3D,sphere,point *,double *);
double garn_pola_cesl(double,double);
void   project_on_sphere(sphere,point,point *);
void   geqn_proj_gotf(point,vect3D,point,point *);
void   cest_reve_fack(c_arc3D,c_arc3D *);
int    qefl_remo_gopd(sphere *,int);
int    hezp_segm_gods(parm,parm,parm,parm,double,double,double,parm *);
double hitf_scal_rikd(vect2D,vect2D);
double rocv_scal_toqc(vect3D,vect3D);
int    segment_circle_inters(parm *,parm,double,parm *);
void   cevg_simp_cafv(int,parm,parm,parm,manif_ro *);
void   guqw_simp_howc(int,parm,parm,parm,
       manif_ro *,int *,int *,int *);
void   lenq_simp_socg(int,int,manif_ro *);
int    melg_spat_nusw(point,double,point,double,
       vect3D,point,double,vect3D,point *);
int    wihz_sphe_vezl(point,double,point,double,double *,
	   point *,double *,double *,double *,double *);
void   vewr_sphe_ruhd(double,double,double,double *,double *);
void   movg_spli_jern(circle3D,point *,int,c_arc3D *);
double qetr_squa_tarj(point,point);
int	   kujw_test_lifn(parm,parm,parm);
int    migz_tole_kums(parm,parm,double);
int    gect_tole_husn(point,point,double);
int    torus_patches(adj_hash,sphere *,int,double,set_arcs *,pt_tor *,int *,int *,int);
int    triangle_circle_inters(parm *,parm,double,c_arc2D *);
void   suzr_find_heqz_wikr(baryc2D,c_arc3D,c_arc3D,
       c_arc3D,point *);
void   qosp_unit_zamk(point,point,point,vect3D *);
void   jatw_unit_hukl(int,manif_ro *);
void   cuwl_unit_pist(parm,parm,vect *);
void   culm_unit_peks(point,point,vect3D *);

void   mejd_allo_dakg(int,int,int,manif_ro *);
void   juqr_allo_nohd(int,int,int,manif_tl *);
void   lagr_allo_goqn(blend_nonself *);
void   foks_allo_vukp(prop_n_curv,ns_curv *);
double ** allocate_mat(int,int);
parm   ** allocate_mat_parm(int,int);
void   homd_allo_tevf(prop_ccurve,c_curve *);
void   gikf_allo_nibl(int,int,int,adj_hash *);

void   fogq_dest_muwf(manif_ro *);
void   zesg_dest_cokd(manif_tl *);
void   newt_dest_lefq(prop_n_curv,ns_curv *);
void   wapl_free_dogc(parm **,int,int);
void   wosn_dest_jomw(prop_ccurve,c_curve *);
void   mewh_dest_selk(int,adj_hash *);
void   vejp_dest_tufq(blend_nonself *);

void   draw_2Dmesh_spec(manif_ro,int);
void   duws_find_vuhk_tosh(manif_ro);
void   zufq_find_nucv_weck(manif_tl);
void   kapj_find_nidr_sotq(manif_tl,point *,int);
void   vuml_find_leks_jupc(manif_tl,point);
void   celv_find_tevg_hejg(manif_tl,rgb_lk *);
void   dufq_find_cepl_qulw(trmsrf *,int);
void   sivh_find_zoql_wefn(parm *,int);
void   digl_find_zuwd_dops(trmsrf *,int *,int);
void   draw_couple_polygons(parm *,int,parm *,int);
void   draw_arc_triangle(parm,double,parm,parm,parm,parm,parm);
void   rogc_find_virw_wufn(pt_tor *,int);
void   loqp_find_mujd_potq(megamanif);
void   dunm_find_roqf_wuln(sphere *,int);
void   qurp_find_dasw_cotr(sphere,manif_tl);
void   qegp_find_wotn_movt(sphere *,manif_tl *,int);
void   gacn_find_qahb_keng(megamanif,PL_curve *,int);
void   gazh_find_cols_hagl(sphere,circle3D *,int);
void   wohs_find_pogf_vezk(sphere *,int,set_arcs *);
void   weqd_find_dahc_kolb(megamanif,PL_curve *,int,point *,int);
void   draw_spheres_and_clouds(sphere *,int,point *,int);
void   cord_find_mecj_nulf(sphere *,manif_tl *,int n,point *,int);
void   rofn_find_sufp_fomn(sphere *,int,set_arcs *,pt_tor *,int);
void   turq_find_tuhn_curv(manif_tl,manif_tl);
void   golc_find_kojr_lewm(trmsrf *,int *,int,c_arc3D *,int);
void   zegd_find_vunh_qubf(pt_tor);
void   draw_treble_arcs(c_arc3D,c_arc3D,c_arc3D);
void   zikt_find_jotz_jewb(trm_sph,trm_sph *);
void   kotg_find_wuhk_kemt(c_curve,c_curve *);
void   jufw_find_lukj_rapg(xform,xform *);
void   jofd_find_mikn_gehj(pt_cnv,pt_cnv *);
void   hepk_find_gict_hubq(circle3D ,circle3D *);
void   zobm_find_wumq_kihf(ns_curv,ns_curv *);
void   getf_find_rogc_todj(point,point *);
void   cunl_find_qedf_rewn(parm,parm *);
void   quvl_find_huts_pagc(manif_tl,manif_tl *);
void   poms_find_resk_lonb(c_arc3D,c_arc3D *);
void   gemq_find_pevg_cikh(PL_curve,PL_curve *);
void   romh_find_cont_qucr(pt_tor,pt_tor *);
void   neqg_find_lodr_bogm(sphere,sphere *);
void   sifm_find_cudw_pafg(blend_nonself,blend_nonself *);
void   goth_find_cofl_futw(c_arc2D,c_arc2D *);
void   guwv_find_dagt_hujw(bd_box3D,bd_box3D *);
 
void   dump_molc_surf(trmsrf *,int ,trmsrf *, int, trmsrf *, int);
