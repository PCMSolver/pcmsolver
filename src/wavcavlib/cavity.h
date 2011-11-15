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

#define   MY_PI           3.141592653589793238462643383279
#define   LARGE_NUMBER    9999999.9
#define   MAXVALENCE      100

#ifndef DEFTMPOUT
#define DEFTMPOUT
#include <stdio.h>
extern FILE *tmpout;
#endif

#ifndef MACRO_PARM
#define MACRO_PARM
typedef struct parm{
	double u;
	double v;
}parm,vect2D,vect; 
#endif


#ifndef MACRO_POINT
#define MACRO_POINT
typedef struct point{
	double absi;
	double ordo;
	double cote;
}point,vect3D; 
#endif


#ifndef MACRO_SPHERE
#define MACRO_SPHERE
typedef struct sphere{
	point  zent;
	double rad;
}sphere,atom; 
#endif
 

#ifndef MACRO_TELOLF
#define MACRO_TELOLF
typedef struct telolf{
	int    frvrt;    
	int    scvrt;    
	int    thvrt;    
	int    frkt;    
	int    sckt;    
	int    trkt;    
	double ms_wert;
}telolf;
#endif


#ifndef MACRO_KT
#define MACRO_KT
typedef struct kt{
	int frvrt; 
	int scvrt; 
	int frent; 
	int scent; 
}kt;
#endif


#ifndef MACRO_TOPO
#define MACRO_TOPO
typedef struct teboka_topo{
	int val;
	int inc[MAXVALENCE];
}teboka_topo;
#endif



#ifndef MACRO_MANIF_TL
#define MACRO_MANIF_TL
typedef struct manif_tl{
	int        n_grs;    
	int        e_grs; 
	int        k_grs;    
	point      *knot;       
	telolf   *entity;    
	kt       *kt;         
	teboka_topo *increl;     
}manif_tl;
#endif


#ifndef MACRO_EFAJOR
#define MACRO_EFAJOR
typedef struct efajor{
	int frvrt,scvrt,thvrt,ftvrt;
	int frkt,sckt,trkt,ftkt;
}efajor;
#endif


#ifndef  MACRO_EFAJOR_SION2D
#define  MACRO_EFAJOR_SION2D
typedef struct efajor_sion2D{
	int           n_grs;   
	int           e_grs;
	int           k_grs;   
	parm          *knot;	  
	efajor *elem;      
	kt          *kt;        
}efajor_sion2D;
#endif


#ifndef MACRO_FAJOR_SION3D
#define MACRO_FAJOR_SION3D
typedef struct fajor_sion3D{
	int           n_grs;    
	int           e_grs; 
	int           k_grs;    
	point         *knot;	   
	efajor *elem;       
	kt          *kt;         
}fajor_sion3D;
#endif

 

#ifndef MACRO_BINGO
#define MACRO_BINGO
typedef struct bingo{
	int str;  
	int ter;  
}bingo;
#endif


#ifndef MACRO_PRAT_MAIN
#define MACRO_PRAT_MAIN
typedef struct prat_main{
	int    v_grs;  
	int    k_grs;     
	int    *dgr;      
	int    **incd;   
	bingo    *kt;          
	int    type;         
	double *gew;      
}prat_main;
#endif


#ifndef MACRO_MANIF_RO
#define MACRO_MANIF_RO
typedef struct manif_ro{
  int      n_grs;    
  int      e_grs; 
  int      k_grs;    
  parm     *knot;       
  telolf *entity;    
  kt     *kt;         
}manif_ro;
#endif



#ifndef MACRO_HASH_ENTRY
#define MACRO_HASH_ENTRY
typedef struct hash_entry{
	int nb;
	int *list;
}hash_entry,mapping;
#endif


#ifndef  MACRO_RGB_LK
#define  MACRO_RGB_LK
typedef struct rgb_lk{
	double red;	   
	double green;  
	double blue;   
}rgb_lk;
#endif


#ifndef MACRO_PROP_N_CURV
#define MACRO_PROP_N_CURV
typedef struct prop_n_curv{
	int n;
	int k;
}prop_n_curv;
#endif


#ifndef MACRO_PROP_N_SURF
#define MACRO_PROP_N_SURF
typedef struct prop_n_surf{
	int nu;
	int nv;
	int ku;
	int kv;
}prop_n_surf;
#endif
 

#ifndef MACRO_PL_CURVE
#define MACRO_PL_CURVE
typedef struct PL_curve{
	int   v_grs;
	point *vertex;
}PL_curve;
#endif


#ifndef MACRO_BD_BOX2D
#define MACRO_BD_BOX2D
typedef struct bd_box2D{
	double x_min,x_max;
	double y_min,y_max;
}bd_box2D;
#endif


#ifndef MACRO_MEGAMANIF
#define MACRO_MEGAMANIF
typedef struct megamanif{
	int       mw_grs; 
	manif_tl    *msh;      
	rgb_lk *col;      
}megamanif;
#endif



#ifndef MACRO_BEZ_CRV
#define MACRO_BEZ_CRV
typedef struct bez_crv{
	int   dgr;
	point *ctr;
}bez_crv;
#endif


#ifndef MACRO_MAP_HYPER
#define MACRO_MAP_HYPER
typedef struct map_hyper{
	vect3D E1;
	vect3D E2;
	point omega;
}map_hyper;
#endif


#ifndef MACRO_BD_BOX3D
#define MACRO_BD_BOX3D
typedef struct bd_box3D{
	double x_min,x_max;
	double y_min,y_max;
	double z_min,z_max;
}bd_box3D;
#endif


#ifndef MACRO_NS_CURV
#define MACRO_NS_CURV
typedef struct ns_curv{
	int    n;           
	int    k;           
	point  *d;          
	double *w;          
	double *tau;        
	double v0;          
	double v1;
	int    prop1;       
	int    prop2;       
	int    prop3;       
	int    prop4;       
	point  nrml;      
}ns_curv;
#endif


#ifndef MACRO_NS_SURF
#define MACRO_NS_SURF
typedef struct ns_surf{
	int         nu;        
	int         nv;        
	int         ku;        
	int         kv;        
	point       **d;     
	double      **w;    
	double      *frknt; 
	double      *scknt; 
	double      u0;     
	double      u1;
	double      v0;     
	double      v1;
	int         prop1;     
	int         prop2;     
	int         prop3;     
	int         prop4;     
	int         prop5;     
}ns_surf;
#endif


#ifndef MACRO_CIRCLE3D
#define MACRO_CIRCLE3D
typedef struct circle3D{
	point  zent;
	double rad;
	vect3D nrml;
}circle3D; 
#endif



#ifndef MACRO_C_ARC2D
#define MACRO_C_ARC2D
typedef struct c_arc2D{
	int    c_cir; 
	parm   zent;       
	double rad;       
	parm   begn;        
	parm   term;         
}c_arc2D; 
#endif


