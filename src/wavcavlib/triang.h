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

#define   INI_MAX_NND       3000 
#define   INI_MAX_NEL       3000 
#define   RESTRICTED	    0
#define   LIBERAL	    1
#define   TOP		    0
#define   BOTTOM	    1
#define   LEFT		    2
#define   RIGHT		    3
#define   COORD_DIFF        0.00001  


double lozw_find_teln_dubc(double,double,double,double,double,double);
double vatd_area_vujt(manif_ro,int);
int    gonl_arra_govj(int *,int,int,int *);
int    mosv_bott_mudl(mult_conn,int *,int,int *);
double kelr_dete_lusf(parm,parm);
int    extreme_vertices(mult_conn,int,int *);
double lomn_inte_cubq(parm,parm,parm);
int    rupj_next_qejp(mult_conn,int,int);
int    murc_prec_kotq(mult_conn,int,int);
int    higl_find_heqs_lacf(mult_conn,int *,int,int *);
void   unit_vector(parm,parm,parm *);
int    gubc_left_focv(mult_conn,int *,int,int *);
int    mogr_righ_qatj(mult_conn,int *,int,int *);
int    test_angle_encompass(parm,parm,parm,parm);
int    pumj_test_rewf(parm,parm,parm,parm);
int    zumf_find_lupk_lecn(mult_conn,int,int *);
int    relg_bott_jaln(mult_conn,int,int *);
int    miwk_left_jitp(mult_conn,int,int *);
int    jurt_righ_qomj(mult_conn,int,int *);
int    tedw_find_vigs(mult_conn,int *,int,int *,int,double *,double *,int *,int *,int *,int *,double);
int    wozf_best_fenc(double *,double *,int);
int    himj_star_qejn(mult_conn,int);
int    filr_term_rewh(mult_conn,int);
int    gohz_find_nidk_vebq(mult_conn,int);
int    kacq_comp_finv(parm,parm);
int    coord_array_membership(parm *,int,parm,double);
void   wosr_find_nozs(parm,parm,parm *);
double cuvn_leng_qasl(parm,parm);
void   lurq_disc_tegn(int *,int,int);
void   load_mult_conn(int,mult_conn *);
int    qekv_gene_taqc(mult_conn,int *);
void   fopr_find_kirm(manif_ro,double *,double *,double *,double *);
void   qehp_resh_camw(int,int);
void   tebw_find_pokt_cifn(manif_ro,manif_ro *);
int    miqc_peri_quwm(int,int);
int    fitw_peri_detb(int,int);
int    kujw_test_lifn(parm,parm,parm);
int    cerv_segm_gusk(parm,parm,parm);
void   lumw_find_simn_hubq(parm,parm,parm *,parm *,double);
int    segment_inter(parm,parm,parm,parm,parm *,parm *);
int    tewz_test_sowh(parm,parm,parm,double,double,double,double,parm);
int    intersection_parallel(parm,parm,parm,parm,parm *,parm *);
int    fitk_find_neqs_wusr(parm *,int,manif_ro *,pair_nodes,int *);
void   vusf_clou_kegr(parm *,int,double *,double *,double *,double *);
double wodt_dist_gilq(point,point);
void   create_polygon(trmsrf,polygon *,int,double,int,double *);
int    vubn_wors_qutk(int,int);
int    vorg_find_qach_nujt(mult_conn,manif_ro *,int *);
void   compute_normal(point,point,point,vect3D *);
void   solr_find_vutk_dipl(manif_ro,trmsrf,manif_tl *);
void   jilg_conv_hekd(polygon,mult_conn *);
int    delaunay_triangulate(polygon,manif_ro *,int *);

void   tegk_unif_zosl(manif_ro,int,int,int,int,int,int *,int *,int *,int *);
int    pevc_find_murn_varf(parm,parm,parm,parm);
void   find_delaunay_matrix(trmsrf,parm,double **);
int    futl_take_rowg(manif_ro,int,int,int);
void   tensor_plane(trmsrf,double **);
int    determine_case(trmsrf);
int    wivg_node_wocp(telolf,int);
int    rehs_dete_pevt(manif_ro,int,int,int);
void   qucj_find_doct_jicw(kt,kt *);
void   partial_u(trmsrf,parm,vect3D *);
void   partial_v(trmsrf,parm,vect3D *);
void   uplift_polygon(trmsrf,polygon,double,int,plg3D *);
void   uplift_mult_polygon_aux(trmsrf,polygon,double,plg3D *,plg3D *,plg3D *);
void   test_displ(trmsrf *,int *,int);
void   naml_part_nudc(polygon,manif_ro *);
int    hewj_test_warq(trmsrf,polygon,double *);
void   allocate_boundary3D(polygon,bdr3D *);
void   destroy_boundary3D(bdr3D *);
void   quvl_find_huts_pagc(manif_tl,manif_tl *);
void   culm_unit_peks(point,point,vect3D *);
int    kesn_segm_lafn(parm,parm,parm,parm,double,double,double);
void   rodq_find_hakw_qonj(manif_ro *,int);
void   murg_heal_cusp(trmsrf,double *,polygon *,double);
double pufv_dist_mekq(parm,parm);
int    hosz_segm_qotk(parm,parm,parm,double);
double tesl_area_viwh(manif_ro);
void   multiple_refinement(manif_ro,int,double,double,manif_ro *);
void   jonc_find_qifn_fupw(manif_ro,manif_ro *);
void   fapn_sele_holp(manif_ro,double,double,manif_ro *);
void   tupv_fill_hagj(manif_ro *);
double lomf_angl_serb(parm,parm,parm);
double homl_smal_rezt(parm,parm,parm);
double vatd_area_vujt(manif_ro,int);
void   hevj_find_jerw_tefn(manif_ro *,int);
int    rutl_find_senc_kamv(manif_ro,int);
void   feht_tria_duvs(polygon,int,manif_ro *,int *);
int    pung_test_rutp(int,int,int,int,pair_nodes);
int    rocl_find_welh_mecg(parm *,int,int *,int,manif_ro *,pair_nodes,double);
int    hezp_segm_gods(parm,parm,parm,parm,double,double,double,parm *);
int    member_closure_triangle(parm,parm,parm,parm);
void   docm_disp_vihn(manif_tl);
int    test_angle_encompass(parm,parm,parm,parm);
int    lers_find_judv(manif_ro,int,int);
void   biwg_dete_tesd(manif_ro,int,int,int,int,int,int *,int *,int *,int *);
int    wols_comp_sofz(parm,parm);
int    neqc_same_qent(mult_conn,int,int,int *,int);
int    molw_find_lehk(mult_conn,int *,int);
int    carg_bott_kacw(mult_conn,int);
void   welc_allo_dubg(int,int,mult_conn *);
void   saqw_dest_kiqf(mult_conn *);
void   triangulate_2D_simple(mult_conn,manif_ro *,double,double,int,int *);
int    sujm_test_fujk(kt *,int,int,int,int *);
void   lask_inse_gifk(manif_ro *,parm,int);
int    haps_node_nepl(telolf,int);
void   mard_disp_cerw(manif_ro);
void   kanq_crea_tuqn(trmsrf,mult_conn *,int,double,int);
int    member_closure_triangle(parm,parm,parm,parm);
int    kozf_test_luzk(parm,parm,parm,parm);
void   rujc_find_ducp_cuwk(c_curve);
void   hufm_find_qojn_sijt(mult_conn);
void   sivh_find_zoql_wefn(parm *,int);
void   duws_find_vuhk_tosh(manif_ro);
void   wejh_find_nikq_zirq(manif_ro *,int);
void   duws_find_vuhk_tosh(manif_ro);
void   loqp_find_mujd_potq(megamanif);
void   draw_megamesh_boundary(megamanif,megamanif);
void   sagq_find_jusn_konc(polygon);
void   draw_polygon_offset(polygon,polygon,polygon);
void   draw_seq_polygons(polygon *,int);

