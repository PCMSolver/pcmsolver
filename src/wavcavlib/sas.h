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

#define INS_SET     0
#define OUT_SET     1
#define INTERF_SET  2
#define Q_TYPE      0
#define T_TYPE      1


#ifndef  MACRO_BLEND_TOR
#define  MACRO_BLEND_TOR
typedef struct blend_tor{
	int trim_idx;   
	int sph_idx1;   
	int sph_idx2;   
	int a_id1;    
	int a_id2;    
}blend_tor;
#endif


#ifndef MACRO_BLEND_CPX
#define MACRO_BLEND_CPX
typedef struct blend_cpx{
	int        bt_grs;
	blend_tor  *BT;         
	hash_entry *HE;         
}blend_cpx;
#endif


#ifndef MACRO_PROP_TOR
#define MACRO_PROP_TOR
typedef struct prop_tor{
	circle3D C1;
	circle3D C2;
	double   D;
	double   r;
	double   R;
	parm     C_pls;
	parm     C_mns;
}prop_tor;
#endif



#ifndef MACRO_PARENT_CIRCLE
#define MACRO_PARENT_CIRCLE
typedef struct parent_circle{
	circle3D supp;    
	int      supp_id; 
	int      sbl; 
}parent_circle;
#endif


#ifndef MACRO_CONC_MODEL
#define MACRO_CONC_MODEL
typedef struct conc_model{
	sphere         S;          
	point          cr0,cr1,cr2;
	manif_tl         msh;        
	point          *G;		   
	vect3D         *nrm;       
	bd_box3D BB;         
}conc_model; 
#endif


#ifndef MACRO_POLYARC
#define MACRO_POLYARC
typedef struct polyarc{
	int            ar_grs; 
	c_arc2D *CA;
}polyarc;
#endif


#ifndef MACRO_NODE_BLEND
#define MACRO_NODE_BLEND
typedef struct node_blend{
	int type_node;    
	int trim_idx;     
	int val;      
	int inc_edge[3];  
}node_blend;
#endif


#ifndef MACRO_KT_BLEND
#define MACRO_KT_BLEND
typedef struct kt_blend{
	int frvrt;
	int scvrt;
	int side_1;
	int side_2;
	int nb_opp;
	int opp[2];
}kt_blend;
#endif


#ifndef MACRO_PRAT_MAIN_BLEND
#define MACRO_PRAT_MAIN_BLEND
typedef struct prat_main_blend{
	int        n_grs;
	int        k_grs;
	node_blend *N;
	kt_blend *E;
}prat_main_blend;
#endif


#ifndef MACRO_INCLUDED_SIDES
#define MACRO_INCLUDED_SIDES
typedef struct included_sides{
	int side_1; 
	int side_2; 
}included_sides;
#endif



#ifndef MACRO_SUPP_SPH_TRI
#define MACRO_SUPP_SPH_TRI
typedef struct supp_sph_tri{
	int sph_idx1;
	int sph_idx2;
	int sph_idx3;
}supp_sph_tri;
#endif


#ifndef MACRO_HALF_INTERS
#define MACRO_HALF_INTERS
typedef struct half_inters{
	int   ts;
	point A0;
	point A1;
}half_inters;
#endif


#ifndef MACRO_MAT_OPERATOR
#define MACRO_MAT_OPERATOR
typedef struct mat_operator{
	double OP[3][3];
}mat_operator;
#endif


int    any_segm_ray(parm,vect2D,parm,parm);
int    hefk_apen_sikf(double,int,sphere *,int,trmsrf *,int,int *,trmsrf *,
	   int,blend_cpx *,set_arcs *,int);
int    cesp_find_qimp_pufv(c_arc3D,point,double);
int    hagf_arcs_belz(adj_hash,double,sphere *,int,set_arcs *,parent_circle *,int);
void   hutw_bbox_zawn(int,double,c_arc3D,bd_box3D *);
void   herg_boun_jesd(megamanif,double *,double *,
       double *,double *,double *,double *);
int    lafc_boun_gusd(bd_box3D,bd_box3D);
void   vupn_chec_wudq(trmsrf *,sphere *,blend_tor);
void   mafj_chec_zogj(c_arc3D,c_arc3D,
	   c_arc3D,c_arc3D);
void   mehz_chec_dijp(int,set_arcs *);
void   check_integrity_blend_complex(trmsrf *,sphere *,blend_cpx,set_arcs *);
int    check_penta_case(trmsrf *,set_arcs *,sphere *,
       int,int,int,int,blend_cpx,c_arc3D,
       point *,point *,point *,c_arc3D *,c_arc3D *,
	   int *,int *,int *,int *);
void   pawt_choo_husn(point,point,c_arc3D *,int,c_arc3D *);
int    jedr_circ_wefj(point,double,circle3D,circle3D,point *);
int    pavz_circ_kuts(circle3D,int,sphere *,int);
void   qumg_circ_tezl(point,point,point,point,double *,point *);
int    qoml_find_venq_fugh(double,sphere *,set_arcs *,int,int,int *,int *,double);
int    refd_comp_lacq(c_arc3D,c_arc3D,
       c_arc3D,trm_sph,c_curve *);
void   sicv_comp_wenr(parm *,int,parm *);
int    connected_by_arc_idx(trmsrf *,blend_cpx,int,int,c_arc3D,double,point *,point *);
void   rifv_find_lips_wecj(mat_operator,mat_operator *);
void   fewd_find_lepz_zaqn(supp_sph_tri,supp_sph_tri *);
void   qavz_find_sowv_tawz(parm *,int,parm,double,parm *);
int    wong_comp_golp(trm_sph,c_arc3D *,int,c_curve *);
int    convert_trim_conc(double,pt_cnv,trmsrf *);
void   pugj_crea_sotm(int,int *,int,trm_sph,c_curve *,int,trmsrf *,int *);
int    tulj_curr_winq(polygon *,int,int *,int);
void   dolj_curv_kacq(trm_sph,c_arc3D,ns_curv *);
void   bikj_deca_vasj(bez_crv,double,point *);
double det_lambda_value(int);
void   direction_spherical_eff(vect3D,mat_operator,vect3D *);
void   kanb_dire_legv(vect3D,mat_operator,vect3D *);
int    directions_int_verif(c_curve,vect2D *);
int    petg_disc_korp(double,c_arc3D);
void   lohs_disc_gojz(c_arc2D,int,parm *);
void   jiwr_disc_qumf(int *,int,int,trmsrf *,int,circle3D,circle3D,set_arcs *,int *);
void   display_BT_id12(int,blend_cpx,set_arcs *);
void   lutc_disp_nulq(c_arc2D);
void   berd_disp_derh(sphere);
double fojn_dist_kerf(point,point,pt_tor,c_arc3D *);
double vuqg_dist_faql(pt_tor,sphere);
double milc_dist_leqk(c_arc3D,c_arc3D);
void   logm_eval_pusn(c_arc2D,double,parm *);
int    numw_edge_jerk(prat_main_blend,int,int);
double rajl_erro_veqd(sphere,c_arc3D);
double mulh_erro_cedm(sphere,point);
double zadt_erro_vehp(sphere,c_arc3D);
void   wolf_eval_murg(trmsrf,double,double,point *);
double jofz_fart_delj(parm *,int,parm);
int    nols_find_lepm_lasv(int,point,double,circle3D,circle3D,pt_tor *,int *);
void   kejh_fill_zogq(sphere *,trmsrf *,int,blend_cpx *);
void   tahq_find_jimp(trmsrf *,trmsrf *,prat_main_blend,double,int *);
void   sedr_fill_luch(int,blend_cpx *);
void   capn_fill_fond(manif_tl *);
int    hicj_fill_zalj(double,sphere *,int,trmsrf *,int,blend_cpx,trmsrf *,supp_sph_tri *,int);
void   fewg_find_duvk(double,sphere,sphere,point *);
void   kotr_find_tujg(manif_ro,int,int,int,pt_cnv,conc_model *);
int    fobw_find_rogs(double,circle3D,circle3D,prop_tor *);
int    find_second_joint(pt_tor,pt_tor,c_arc3D,point,point,point *,point *);
void   volq_find_luvc_jurq(trmsrf *,int);
int    qonk_form_sevt(int,trm_sph,int,c_curve *,int,trmsrf *,int,int *,int);
int    gap_trim_surf(adj_hash,trmsrf *,int,
       set_arcs *,sphere *,int,blend_cpx *,int);
void   dels_geod_tuzd(point,double,point,point,c_arc3D *);
void   cojk_grad_negl(parm *,int,parm,int,manif_ro *);
int    luhb_inco_zokt(parent_circle *,adj_hash,double,sphere *,
	   int,set_arcs *,int,trmsrf *,blend_cpx *,int);
int    ticj_inte_ludj(conc_model,conc_model);
int    qiwj_inte_qivs(parm *,double,circle2D *,circle2D *);
void   qoml_inve_fasd(vect3D,mat_operator,vect3D *);
void   raqn_inve_kotj(int,set_arcs *);
void   involved_entities_BT(trmsrf *,adj_hash,set_arcs *,
       blend_cpx,int,int *,int *,c_arc3D *,int *,int *,
       int *,int *,int,int,int);
int    wucl_find_henl_lorh(trmsrf);
double biqf_leng_lozs(c_arc2D);
int    mefj_list_hepq(parent_circle *,adj_hash,int,int,sphere *,int,double,set_arcs *,int *,int *,int *,int *);
int    hect_list_mudc(sphere *,int,adj_hash,set_arcs *,trmsrf *,int *,int);
int    list_nond_conn(sphere *,trmsrf *,blend_cpx,set_arcs *,
       int,int,c_arc3D *,int *);
int    gazq_loca_joth(sphere *,trmsrf *,int,blend_cpx,
       c_arc3D *,c_arc3D *,c_arc3D *,supp_sph_tri *);
void   kihv_find_firn_qivw(double,double,mat_operator *);
void   surq_find_tejq_kedm(double,double,mat_operator *);
void   fisd_find_rumj_veft(double **,double *,double *);
double fewq_find_qukn_qows(double *,int);
void   kecv_mesh_honc(int,int,int,int,manif_ro *,int *);
void   puwj_midp_curq(c_arc3D,point *);
void   zikf_midp_kusd(c_arc3D,mat_operator,
       mat_operator,point *);
void   nuvp_modi_qetl(trmsrf *);
void   kong_most_bejc(parm *,int,int *,int *,int *);
void   zisq_find_cowk_zevq(int,sphere,sphere,point,circle3D,circle3D,blend_nonself *);
int    cezv_obso_nofm(int,circle3D,parent_circle *,int,int *);
void   suwf_orga_cezm(double,c_arc3D *,c_arc3D *,c_arc3D *,c_arc3D *);
void   orthogonal_center(point,point,point,point,double,
       double,double,double,point *);
void   cerw_orga_sevd(c_arc3D *,c_arc3D *,
       c_arc3D *,c_arc3D *,c_arc3D *);
void   gojw_pair_nokv(int,set_arcs *,parent_circle *,int);
int    vejg_para_rilm(double,sphere,sphere,circle3D *,circle3D *);
void   hevb_para_tucp(double,sphere,sphere,circle3D *,circle3D *);
void   jotc_find_jewl_zorv(int,prop_ccurve *);
void   lasn_find_vupn_pehr(int,prop_ccurve *);
void   mesl_punc_guvf(circle3D,double,point *);
void   sovt_plan_tocz(circle3D,point,point,point,point *);
void   jets_find_cifl_kitq(c_curve,polyarc *);
int    vahg_find_purt_tumh(parm,c_arc2D *,int);
int    gedh_poss_jucp(sphere *,double,int,int,set_arcs *,
       int *,int *,int,int,pt_tor *,int *,int *,int);
double rand_unit();
int    ray_inclusivity(parm *,int,polygon,vect2D);
int    relevant_arc_tor(double,sphere *,int,c_arc3D,
       int,adj_hash,prop_tor *,int);
int    zosd_remo_qohj(atom *,int);
int    rem_penta_patches(double,trmsrf *,int,set_arcs *,
	   sphere *,int,blend_cpx *,int);
double mijt_scal_letp(int,double *,double);
void   serp_scal_belt(megamanif *,double,double,double,double,double,double);
void   toqn_shif_cifm(megamanif *,double,double,double,double,double,double);
int    secray_inclusivity(polygon,polygon);
void   sofl_segm_salc(trmsrf,c_curve,point *);
void   fiwh_simm_tucs(int,manif_ro *,int *,int *,int *);
void   wocf_simp_fird(parm *,int,parm,int,manif_ro *);
int    sapn_some_judq(int,trmsrf *,int,int *,parent_circle *,int *,
	   adj_hash,sphere *,int,double,set_arcs *,trmsrf *,int,
	   int,int,int,int *,int *,double,int *);
void   homg_sphe_vemp(double,trmsrf *,int);
int    sahf_spli_gehq(c_arc3D,point *,int,c_arc3D *);
void   dept_spli_piwg(circle3D,double,double,
       point *,int,c_arc3D *);
void   rezc_spli_qizk(circle3D,mat_operator,
       mat_operator,point *,int,c_arc3D *);
int    pubd_spre_pedc(trmsrf *,trmsrf *,prat_main_blend,
       int *,int *,double);
int    hokq_stri_jipg(sphere,point,double);
void   zevt_stru_copb(adj_hash H,atom *,int,set_arcs *,int *,
	   trmsrf *,int,trmsrf *,blend_cpx,int *);
int    fuqc_test_kuml(c_arc3D,c_arc3D,
       double,double,double,double);
int    touching_BT(trmsrf *,blend_cpx,int,int,double);
void   lekr_find_sujr_niqc(baryc2D,c_curve,parm *);
int    paqf_trim_cist(double,c_arc3D,
       c_arc3D,c_arc3D,trm_sph *,double);
void   sonv_trim_gopl(sphere *,int,trm_sph *,adj_hash);
void   update_merge(trmsrf *,set_arcs *,int,blend_cpx *,parent_circle *);
void   update_set_arcs(set_arcs *);
void   legw_upda_tevg(int,blend_cpx *,int *,int *,int);
void   wasd_veri_vurc(sphere *,trmsrf *,
       trmsrf *,int,blend_cpx,supp_sph_tri *);
void   vowc_veri_cotq(trmsrf *,trmsrf *,prat_main_blend);

double distance_pat_tor_point(pt_tor,point,point *);
double distance_pat_toroidals(pt_tor,pt_tor,point *,point *);
double jerw_dist_zojc(pt_tor,point,point,c_arc3D *);
double hosf_dist_jocw(pt_tor,c_arc3D);
double distance_pat_tor_arc(pt_tor,c_arc3D,point *);

int    peqj_test_vesf(circle3D,sphere);
void   test_convertion_concave(double,trmsrf *,int);
int    test_coparent_arcs(c_arc3D,c_arc3D);
int    farw_test_fijh(plane *,int,int,point);
int    lupt_test_fung(c_arc3D,circle3D,double,double,double);
int    tevb_test_vefp(c_arc3D,circle3D);

void   qatn_allo_mobc(int,blend_cpx *);
void   hegp_allo_qogc(trmsrf *,int,int);
void   culj_allo_pudn(trmsrf *,int);
void   gilp_allo_temc(trmsrf *,int);

void   bulf_dest_tacs(int,blend_cpx *);
void   zaqw_dest_jelq(trmsrf *,int,int,int);
void   novw_dest_wojm(trmsrf *,int,int);
void   luqw_dest_horb(trmsrf *,int,int);

void   fulc_find_lefp_norh(sphere *,int,set_arcs *);
void   tokf_find_javr_pinm(int,int,sphere *,set_arcs *,adj_hash);
void   fiwd_find_qokz_sutv(int,sphere *,set_arcs *,adj_hash);
void   kahr_find_qewr_cows(sphere *,int,pt_tor *,int);
void   levr_find_hulb_hezm(trmsrf *,int,trmsrf *,
	   int,trmsrf *,int);
void   nuvm_find_rewl_pufq(parm *,int);
void   mict_find_juwq_hemq(parm *,int,parm);
void   purm_find_guzf_tecf(point,point,point,double,double,double);
void   lokt_shif_zuhf(PL_curve *,int,double,double,
       double,double,double,double);

