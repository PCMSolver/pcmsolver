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

#define   WEIGHT_LIGHT    0.1      
#define   WEIGHT_HEAVY    2.0      
#define   WEIGHT_SHARP    7.0      
#define   SML_DEPTH       1.5
#define   LRG_DEPTH       6.0
#define   COLOR_LIGHT     0
#define   COLOR_HEAVY     1
#define   COLOR_SHARP     2
#define   PROGRESS        1
#define   STAGNATION      0


#ifndef MACRO_FUND_QUAD
#define MACRO_FUND_QUAD
typedef struct fund_quad{
	double A[4][4];
	double b[3];
	double c;
}fund_quad;
#endif


#ifndef MACRO_PL_EQ
#define MACRO_PL_EQ
typedef struct pl_eq{
	double a;
	double b;
	double c;
	double d;
}pl_eq;
#endif


#ifndef MACRO_SPAN_HZ
#define MACRO_SPAN_HZ
typedef struct span_hz{
	int wrz;              
	int v_grs;       
	int *v_value;         
	int *nb_children;      
	int **chd;        
	int *prt;           
	int dp;             
	int *niv_length;     
	int **niv;           
	int *pos;         
}span_hz; 
#endif


void   sejm_allo_bump(int,span_hz *);
void   mokq_chec_nukb(manif_tl);
void   pajc_coll_qetv(manif_tl *,int,point,int *,int *);
void   roph_coll_fact(manif_tl *,int,int,point);
void   maqb_comp_sogf(manif_tl,manif_tl *);
int    lufj_cons_retw(manif_tl,int,int,point);
void   vegs_disc_lepq(manif_tl *,int *,int);
void   topr_disc_nufr(manif_tl *,int,int *);
void   fomd_deci_todj(manif_tl,manif_tl *,double,double,int,int);
void   qemt_dire_lajk(manif_tl,int,fund_quad *);
double wohn_dire_hurg(manif_tl,int,int,point *);
void   sujt_deal_wozc(int,span_hz *);
void   vocm_disp_qosh(span_hz);
int    gubs_feas_dujz(manif_tl,int,double,point);
void   detl_find_faqs_cakq(fajor_sion3D *,int,int *);
void   ceqj_fund_jolg(manif_tl,int,fund_quad *,fund_quad *);
void   dopr_fuse_tenv(manif_tl *,int,int,int,int,int,int,int,int,int *);
int    sokm_loca_vazn(manif_tl,int);
int    cugt_loca_pogt(manif_tl,int);
int    darj_loca_wejf(manif_tl,int);
int    nuwq_mini_solr(fund_quad,point *);
double lirn_find_qokc_tecj(prat_main,int,span_hz *);
double selr_opti_zoqp(manif_tl,int,int,fund_quad,fund_quad,point *);
void   kodr_plan_pifn(point,point,point,pl_eq *);
void   vowg_plan_vowt(point,point,point,pl_eq *);
int    balt_posi_wurp(point,point,pl_eq);
void   lomr_quad_gudw(span_hz,int,manif_tl,fajor_sion3D *,int,int);
void   decs_quad_sibl(span_hz,int,manif_tl *,fajor_sion3D *,int,int);
double satf_quad_quzl(fund_quad,point);
void   jasg_rear_cupg(int,double *,int *);
void   ciql_remo_cekq(fajor_sion3D *,double,int);
void   tups_remo_veqk(fajor_sion3D *);
void   fuwr_repl_wejf(manif_tl *,int,int);
void   conl_sing_cidp(manif_tl,int,fund_quad *);
void   gedt_telo_fesv(span_hz *,manif_tl *,fajor_sion3D *,int,int);
void   qetk_trim_zogf(span_hz *,int);
void   ruqn_upda_ducz(manif_tl *,int *,int *,
       int *,int,int *,int,int,double *,fund_quad *);


