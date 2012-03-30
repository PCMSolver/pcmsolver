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

#define   CHANGE_ID       123.45679
#define   ADJACENT_PAIR   0
#define   DISJOINT_PAIR   1
#define   SUCCESS         1
#define   FAILURE         0
#define   TYP_COAR        0
#define   TYP_FIN         1


void levp_adap_nirz(double *, int, int, int *);
double potb_anis_dopg(point, point, point);
double valp_area_qelk(point, point, point);
void begin_cell();
void gewn_chec_qamk(manif_tl);
void hetc_conv_lujn(manif_tl, fajor_sion3D *, int, int);
void convert_routine(manif_tl);
int dosc_coor_licf(point *, int, point, double, int *);
void gohl_disc_gunp(manif_tl *, int);
void zudp_disc_duwc(manif_tl *, int, int *);
void nusp_disc_jidk(double, double, manif_tl *, rgb_lk *, sphere *, int);
void wihp_disc_sogj(double, double, manif_tl *, int);
void tugw_extr_vuzf(manif_tl, manif_tl *, int *, int, int *, int *, int, int);
void hanw_fill_keph(prat_main *, int);
void hobf_fill_cewr(manif_tl *);
int cefg_find_fovw(manif_tl, int, int);
void goft_find_doqk(manif_tl, double *, double *, double *, double *, double *, double *);
void teqr_fill_regm(fajor_sion3D *, int);
void cogv_fill_zicd(manif_tl *, int);
void qosr_fill_fedt(manif_tl *);
int cepj_inci_jeqt(manif_tl, int, int *);
double gorh_inte_mesf(point, point, point);
void pebg_load_dogf(int, int, fajor_sion3D *);
void gopd_load_puqs(int, int, fajor_sion3D *);
void taqr_mesh_tuqf(manif_tl, prat_main *, int *, int *);
int relh_numb_nulw(int);
int cuzh_next_kadm(fajor_sion3D, int, int);
int zoqw_numb_lohk(int);
void nahd_find_vukd_gaph(manif_tl *, double, double, double, double, double, double);
int vefd_test_wacj(parm, parm, parm, parm);
void tujh_orie_novd(manif_tl *);
int kuts_next_kolv(telolf, int);
void furk_prep_huqs(manif_tl *, int);
void dopj_remo_qosd(manif_tl *);
int husr_proj_qukw(point, manif_tl);
double rocq_stat_todr(double, manif_tl, int, int, int *, int *);
void gekj_tran_rift(manif_tl, manif_tl *);
void test_ord();
void poml_trid_fert(double **, double *, double *, int);
double qetr_squa_tarj(point, point);
void hegc_orie_laqn(fajor_sion3D *, int);
int cojs_shar_sejm(fajor_sion3D, int, int, int *);
void vetk_stor_lunr(int, fajor_sion3D);
void pesw_stor_neqh(int, fajor_sion3D);
void gows_veri_jiqw(fajor_sion3D);
int qulw_vois_nuhk(manif_tl, int, int, int *);
void jacm_chec_lejt(manif_tl);
void wipz_inci_gemk(manif_tl, hash_entry *);
int vist_find_jetl(manif_tl, int, int *, hash_entry *);
void wacl_coll_hijt(manif_tl *, int *, int, int *, int *);
int fahg_edge_lujd(double, double, manif_tl, int *, int, int *, int *, int);

void lawn_dest_jukt(manif_tl *);
void conw_dest_vojk(int, prat_main *);

void vogn_allo_cusr(int, int, int, prat_main *);
void rudk_allo_tamq(int, int, int, manif_tl *);

void qiwc_find_weqt_cuqh(point *, int);
void lenq_find_devz_wasj(manif_tl, int, int);
void meqr_find_jolt_sobr(manif_tl, point *, point *, int);
void mejn_find_merj_cawt(manif_tl *, int *, int);
void fevg_find_ripd_qeft(fajor_sion3D);
void rebc_find_jetn_jadh(manif_tl, prat_main);
