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

#define ON_ALPHA    0
#define ON_BETA     1
#define ON_GAMMA    2
#define ON_DELTA    3
#define START_NODE  0
#define TERMI_NODE  1
 


#ifndef MACRO_SP_ND_M
#define MACRO_SP_ND_M
typedef struct sp_nd_m{
	int idx;       
	int surf_set;  
	int c_grs;
	int *corner;   
	int val;    
	int *inc_edge_m;  
}sp_nd_m;
#endif


#ifndef MACRO_KT_M
#define MACRO_KT_M
typedef struct kt_m{
	int pt1;       
	int pt2;       
	int corner1[2];
	int corner2[2];
}kt_m;
#endif


#ifndef MACRO_PRAT_MAIN_M
#define MACRO_PRAT_MAIN_M
typedef struct prat_main_m{
	int    n_grs;
	int    k_grs;
	sp_nd_m *N;      
	kt_m *E;
}prat_main_m;
#endif


#ifndef MACRO_SUPP_SURF
#define MACRO_SUPP_SURF
typedef struct supp_surf{
	int s_type;
	int s_id;
}supp_surf;
#endif


#ifndef MACRO_KT_CORN
#define MACRO_KT_CORN
typedef struct kt_corn{
	int            pos; 
	int            n1;  
	int            n2;  
	bd_box3D BB;
	point          mid_node;
}kt_corn;
#endif


#ifndef MACRO_MSH_CORN
#define MACRO_MSH_CORN
typedef struct msh_corn{
	int		    h_grs; 
	int		    k_grs; 
	hash_entry  corn_ext; 
	hash_entry  *corn_int;
	kt_corn   *E;       
}msh_corn;
#endif


#ifndef MACRO_INCIDENT_GNODE
#define MACRO_INCIDENT_GNODE
typedef struct incident_gnode{
	int msh_id; 
	int kt_id;  
}incident_gnode;
#endif


#ifndef MACRO_CONN_SEQ
#define MACRO_CONN_SEQ
typedef struct conn_seq{
	int        nb_comp;
	hash_entry *seq;   
}conn_seq;
#endif


#ifndef MACRO_FAM_NODES
#define MACRO_FAM_NODES
typedef struct fam_nodes{
	int   n_grs;
	int   *msh_id;
	int   *knot_id;
	point X;
}fam_nodes;
#endif


#ifndef MACRO_PROP_DISCR
#define MACRO_PROP_DISCR
typedef struct prop_discr{
	int    dp;
	double id_len;
	double fehl;
}prop_discr;
#endif



#ifndef MACRO_SIDE_FLAG
#define MACRO_SIDE_FLAG
typedef struct side_flag{
	int  aflg;
	int  bflg;
	int  gflg;
	int  dflg;
}side_flag;
#endif



#ifndef MACRO_EDGE_MC
#define MACRO_EDGE_MC
typedef struct edge_mc{
	int n_str;  
	int n_ter;
	int id_cpt;
	int pos;   
}edge_mc;
#endif



#ifndef MACRO_ADD_MC
#define MACRO_ADD_MC
typedef struct add_mc{
	int     nb_edge_mc;
	edge_mc *EMC;
}add_mc;
#endif



#ifndef MACRO_ARC_SUPP
#define MACRO_ARC_SUPP
typedef struct arc_supp{
	int n_arcs;
	int *q_sa;
}arc_supp;
#endif



#ifndef MACRO_SUPPORTING
#define MACRO_SUPPORTING
typedef struct supporting{
	int nin;
	int *supp_ext;
	int **supp_int;
}supporting;
#endif


int    adaptive_rectif(atom *,int,blend_cpx,set_arcs *,arc_supp *,
	   trmsrf *,int,int *,trmsrf *,int *,int,int *,int *);
void   allocate_add_mc(int,add_mc *);
void   pucb_allo_pakq(msh_corn,msh_corn *);
void   jizw_allo_vipg(trmsrf *,int,int,int,int,prat_main_m *);
void   arc_disc_param(blend_cpx,trmsrf *,int,set_arcs,arc_supp *);
int    tanv_asso_qupr(int,int,megamanif,int *,msh_corn *,int *,int *);
int    welq_asso_benv(int,int,megamanif,int *,msh_corn *,point **,int *,int *,int **);
void   luvm_blen_foms(atom *,int,int *,trmsrf *,int,trmsrf *,
       int,trmsrf *,int,set_arcs *,blend_cpx,hash_entry *,
	   int,double,int *);
void   hojp_boun_hapw(manif_tl,conn_seq *,int);
void   letd_boun_noth(pt_tor,double,bd_box3D *);
void   nojp_boun_zovq(pt_cnv,double,bd_box3D *);
void   hajk_chec_wizn(manif_tl);
int    check_mesh_integrity_forc(manif_tl);
int    fizt_chec_fuwk(mult_conn,double,parm *);
void   fenq_find_zoft_daqr(msh_corn,msh_corn *);
void   kulf_dest_lajs(int,int,int,int,prat_main_m *);
void   destroy_add_mc(add_mc *);
int    logj_dete_nehl(manif_tl);
void   det_graph_bl_simple(sphere *,trmsrf *,int,trmsrf *,
	   int,blend_cpx,prat_main_blend *G,included_sides *inc);
void   lutn_dest_mogt(msh_corn *);
void   tuwq_disc_nafp(double,msh_corn *,point **,megamanif *,int *,double,double,int *);
void   muhs_disc_kesq(sphere *,trmsrf,trmsrf *,
       blend_cpx,int,set_arcs,int *,int *,int **,int *,arc_supp);
void   sazk_disc_tujc(prop_discr,prop_discr,atom *,int,int *,trmsrf *,
	   int,trmsrf *,int,trmsrf *,int,blend_cpx,set_arcs *,megamanif *,
	   megamanif *,supp_surf *,supp_surf *,prat_main_m *,prat_main_m *,int,msh_corn *,msh_corn *,
	   point **,point **,int *);
void   hift_disp_tevr(msh_corn);
void   display_side(int);
void   vijm_disp_rofn(int);
double qehs_dist_zugt(int,int,msh_corn,msh_corn,manif_tl,manif_tl);
double wujt_dist_kaqt(side_flag,pt_cnv,
       c_arc3D,int *,int *);
double jufq_dist_dunt(side_flag,pt_tor,
       c_arc3D,int *,int *);
int    edges_coarse_ext(mult_conn,add_mc,int *,int *);
int    edges_self_inters(mult_conn,add_mc,int *);
void   wegf_find_hord_luvk(manif_tl);
void   somt_find_qoph_hegd(manif_tl,rgb_lk *,trmsrf *,
       trmsrf *,trmsrf *,supp_surf *);
void   deqt_find_vegc_cobn(manif_tl);
void   dulw_find_tigv_dubr(manif_tl);
int    qunf_expa_foqg(manif_tl *,int *,hash_entry *);
void   wovn_fill_nelr(double,manif_tl,point *,msh_corn *);
void   cogv_fill_zicd(manif_tl *,int);
void   tuwh_fill_zupf(msh_corn,int,int,int,double,manif_tl,msh_corn *,int *,int *); 
void   lurp_find_kafg_delc(double,int,manif_tl *,msh_corn *,int,point **,int *);
void   hucz_fuse_qorw(int,int,manif_tl *,int *,int *);
int    jish_inci_gelh(sphere *,trmsrf *,
       blend_cpx,int,int *,int *,double);
int    id_blend_for_coarse(trmsrf,supporting,
       mult_conn,edge_mc,set_arcs,arc_supp);
void   nohv_find_zatg_meqt(pt_tor,int,int,c_arc3D *);
void   kumn_inci_fuzr(manif_tl,hash_entry *);
int    wuvp_inci_tejg(megamanif,int *,msh_corn *,point **,incident_gnode **);
void   initialize_disc_len(int,double,trmsrf *,prat_main_blend,
	   int *,included_sides *inc,int *);
void   zekc_merg_gelw(megamanif *,msh_corn *,point **,supp_surf *,supp_surf *,rgb_lk *,int,manif_tl *,int *);
int    mult_conn_simple(trmsrf,int *,
       int **,mult_conn *,add_mc *,int);
void   gick_mult_jufs(int *,pt_cnv,mult_conn *);
int    juzv_mult_gehv(trmsrf,int *,int **,mult_conn *,add_mc *,int *,msh_corn *,point *);
int    kuts_next_kolv(telolf,int);
void   cisv_find_kuns_fenj(trmsrf,mult_conn,
       manif_ro *,double,int,int,int,int,int *);
int    nosp_find_nozm_woqr(manif_ro,int,int);
int    nevl_oppo_cikj(prat_main_blend,int,int *);
void   kecq_orie_watk(manif_tl *,int);
void   somn_popu_wolm(manif_tl,manif_tl *);
void   qazf_popu_fawj(manif_tl,rgb_lk *,trmsrf *,
       trmsrf *,trmsrf *,supp_surf *,manif_tl *,
	   rgb_lk *,sphere *);
void   hisl_popu_nuvw(manif_tl,manif_tl *);
void   pows_popu_peql(manif_tl,manif_tl *);
int    digv_pour_newl(trmsrf,int *,int **);
void   najt_prop_pics(megamanif,int *,int *,int *);
void   qedc_smoo_loct(manif_ro *,int,int*);
int    wehc_sphe_bufp(trmsrf,sphere,mult_conn,
       manif_ro *,double,int,int,int,int,int*);
void   somn_take_dogm(int,msh_corn,int *,int *);
int    kotc_test_kolj(manif_tl,int,int);
int    salr_test_jofl(trmsrf,sphere,double,double *);
int    lobj_thin_fizn(manif_tl,double);
void   munk_upda_vugz(int,int,int *,int,msh_corn *);
void   peql_mult_tenq(trmsrf surf,sphere S,
       manif_ro mshin,int max,double accuracy,manif_ro *mshout,int pat,
       int mx_nnd,int mx_nel,int mx_ned);
void   cetm_find_cojq_vurg(trmsrf *,trmsrf *,
       bd_box3D *,prat_main_blend *,int *,included_sides *,
       double,double,int);

