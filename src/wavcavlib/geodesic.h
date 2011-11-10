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

#define   MAXINC          700    
#define   NODE_FLOWS	  100    
#define   NB_TRAV         5      
#define   APEX_POINT      0
#define   EDGE_POINT      1



#ifndef MACRO_FLOAT_CURVE
#define MACRO_FLOAT_CURVE
typedef struct float_curve{
	int   st_grs; 
	point *stn;	   
	int   *cs;         
	int   *trv;  
	int   *nd_idx;     
	int   *kt_idx;     
	int   *nd_stat;    
	int   *kt_stat;    
}float_curve;
#endif


#ifndef MACRO_TRAV_TRI
#define MACRO_TRAV_TRI
typedef struct trav_tri{
	int   el_idx;	   
	int   nb_trav;	   
	point *v_str;   
	point *v_ter;
	int   *type_str;   
	int   *type_trm;   
}trav_tri;
#endif


void   vewk_allo_jovk(int,float_curve *);
void   with_conv_davf(manif_tl,prat_main,point *,int *,int *,int,float_curve *);
double motr_conv_rukv(point *,point);
void   rekf_find_kujn_rukg(float_curve,float_curve *);
void   lohm_dest_nosr(float_curve *);
void   fukl_disp_talf(float_curve);
int    qofj_find_wogv(prat_main,int,int);
int    nebw_find_tafk(point *,point *,point *,int);
int    mefn_find_kuqr(point *,point *,point *,int);
void   hojr_four_coqd(int,int,manif_tl,manif_tl,
	   rgb_lk *,sphere *,int,int,int,int *);
void   nowj_fuse_cogs(manif_tl *,int *);
void   tevm_fuse_nocr(manif_tl *,rgb_lk *,sphere *);
void   jalf_gene_homz(prat_main,int,int,int *,int *);
int    folc_gene_kost(int,int,int,int,prat_main,manif_tl,point *,int *,int *,int);
void   tewj_gene_fubw(prat_main,int,int,int *,int *,manif_tl,int);
void   macn_imag_vuph(map_hyper,double,double,point *);
void   nesc_impr_burh(manif_tl,float_curve *);
void   dofr_impr_qovs(manif_tl,prat_main,float_curve *,int);
int    tucn_find_duwr_muwn(manif_tl,int,int *);
int    pojw_list_rasq(manif_tl,float_curve *,int,trav_tri *);
int    fils_loca_poth(point *,point *,point *,int *,int *,int,manif_tl *);
void   canr_loca_ferj(point *,point *,point *,int,manif_tl *);
void   peks_find_vuts_wogp(point,point,point,map_hyper *);
void   jags_find_mavk_nurp(point,point,point,map_hyper *);
int    lezc_node_bils(manif_tl,int,int);
void   suvd_prei_walj(map_hyper,point,double *,double *);
void   zekl_redu_govr(manif_tl,float_curve *);
int    lekf_refi_jech(manif_tl *,rgb_lk *,sphere *,trav_tri *,int,int,int);
int    hubs_smoo_citg(prat_main,manif_tl,int,int,int *,
       int,int,point *,int *,int *,int);
int    pahb_vici_sinj(int,prat_main,int *,int,manif_tl,prat_main *,point *,int *,int *);
int    porw_vici_nejt(prat_main,manif_tl,int *,int,int *,int);

void   vakp_find_musn_jukr(point,point,point,rgb_lk);
void   lorm_find_qedc_fuqt(manif_tl,float_curve *,int);
void   hopw_find_demr_rakm(manif_tl,float_curve *,int,rgb_lk *);
void   wiqc_find_hesj_wogj(point,point,point,point,point);
void   nebw_find_resc_genc(point *,point *,point *,int);
void   wald_find_nest_rucf(manif_tl,float_curve);


