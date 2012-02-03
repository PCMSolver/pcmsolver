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

#define  NORTH_HEMI  0
#define  SOUTH_HEMI  1


#ifndef  MACRO_QUADRANGULATION
#define  MACRO_QUADRANGULATION
typedef struct quadrangulation{
	int			  n_grs;   
	int			  e_grs;
	int			  k_grs;   
	parm          *knot;	  
	double        *zt;      
	int           *flag;      
	int           *type;	  
	efajor *elem;      
	kt          *kt;        
}quadrangulation;
#endif


void  sotp_atom_kegt(atom *,int);
void  rahv_atom_zelf(trmsrf,quadrangulation *);
void  atom_iso_once(FILE *,atom,int);
void  kucw_atom_nocd(int,atom,trmsrf *);
int   qumw_atom_picm(atom *,double,trmsrf *,int,
      trmsrf *,int,trmsrf *,int,set_arcs *,
	  blend_cpx,double,int,hash_entry *,int*,int,int);
int   gikl_atom_vumh(double,c_arc3D *,
      int,trmsrf *);
int   qidk_conn_bowz(c_arc3D *,
      int,hash_entry *,double);
int   varm_dete_qumz(atom *,int,double,int *);
void  caqw_disp_pegr(prat_main_blend);
void  ditq_find_jekr_gevw(c_arc3D,c_arc3D,c_arc3D);
void  itl_coons_patch(trmsrf,
      double,double,parm *);
void  corm_fill_ruqt(quadrangulation *);
void  gojw_quad_wuln(int,trmsrf,
      quadrangulation,int,c_curve *);
void  trimmed_probe_hole(c_arc3D *,int,
      trm_sph,trmsrf *);
int   lern_unat_moln(trmsrf *,prat_main_blend,
      int,c_arc3D *,int *);




