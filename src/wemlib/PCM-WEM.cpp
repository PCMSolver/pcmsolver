#include <iostream>
#include <fstream>
#include <cmath>

#include "PCM-WEM.h"
// Beautiful...
extern "C"{
#include "WEM.h"
#include "read_points.h"
#include "vector2.h"
#include "interpolate.h"
#include "topology.h"
#include "kern.h"
#include "compression.h"
#include "postproc.h"
#include "WEMRHS.h"
#include "WEMPCG.h"
#include "WEMPGMRES.h"
#include "dwt.h"
#include "cubature.h"
#include "gauss_square.h"
#include "constants.h"
}


// Greens function wrappers for Helmut's C-code

static GreensFunctionInterface *gf = NULL;

static double SingleLayer(vector3 x, vector3 y){
  
  Vector3d vx(x.x, x.y, x.z);
  Vector3d vy(y.x, y.y, y.z);
  
  return gf->getSingleLayer(vx,vy);
}

static double DoubleLayer(vector3 x, vector3 y, vector3 n_y){
  
  Vector3d vx(x.x, x.y, x.z);
  Vector3d vy(y.x, y.y, y.z);
  Vector3d vn_y(n_y.x, n_y.y, n_y.z);
  
  return gf->getDoubleLayer(vx,vy,vn_y);
}

WaveletBEM::WaveletBEM(){

  P_ = NULL;
  F_ = NULL;
  T_ = NULL;
  E_ = NULL;
  W_ = NULL;
  systemMatricesInitialized_ = false;

  eps_ = 1e-10;
  quadratureLevel_ = 1;

  nPoints_ = 0;
}

WaveletBEM::~WaveletBEM(){
  if(P_ != NULL) free(P_);
  if(F_ != NULL) free_patchlist(&F_,nf_);
  if(T_ != NULL) free_interpolate(&T_,p_,M_);
  if(E_ != NULL) free_elementlist(&E_,p_,M_);
  if(W_ != NULL) free_waveletlist(&W_,p_,M_);

  if(systemMatricesInitialized_){
    free_sparse2(&S_i_);
    free_sparse2(&S_e_);
  }
}

BEMInterface *WaveletBEM::newBEMInterface(){
  return new WaveletBEM();
}

bool WaveletBEM::readCavity(std::string &filename){

  vector3 ***U;

  int res = read_points1(&U, &p_, &M_, const_cast<char *>(filename.c_str()));
  if(res != 0){
    // Or better, throw something
    //std::cerr << "Error opening cavity file \"" << filename << "\"" << std::endl;
    return true;
  }

  nf_ = p_*(1<<M_)*(1<<M_);
  /* Topologie bestimmen */
  init_interpolate(&T_,U,p_,M_);
  np_ = gennet(&P_,&F_,U,p_,M_);

  free_points(&U,p_,M_);

  // Calculate volume and area here, and also a list of quadrature
  // points and normals

  return false;
}


void WaveletBEM::getCavityDefs(std::vector<Vector3d> &points,
			       std::vector<Vector3d> &normals,
			       std::vector<double> &areas){  

  vector2 s,t;
  vector3 point;
  vector3 normal;
  unsigned int n = 1 << M_;
  double h = 1./n;
  cubature *Q;
  unsigned int g = quadratureLevel_;
  init_Gauss_Square(&Q,g+1);

  for (unsigned int i1=0; i1<p_; i1++){
    s.y = 0;
    for (unsigned int i2=0; i2<n; i2++){
      s.x = 0;
      for (unsigned int i3=0; i3<n; i3++){
	
	for (unsigned int k=0; k<Q[g].nop; k++){
	  t = vector2_add(s,vector2_Smul(h,Q[g].xi[k]));

	  point = Chi(t,T_[i1],M_);
	  Vector3d po(point.x, point.y, point.z);
	  points.push_back(po);	
  
	  normal = n_Chi(t,T_[i1],M_);
	  Vector3d no(normal.x, normal.y, normal.z);
	  normals.push_back(no);
	  	  
	  double area = h*h*Q[g].w[k]*vector3_norm(n_Chi(t,T_[i1],M_));
	  areas.push_back(area);
	  }
	
	s.x += h;
      }
      s.y += h;
    }
  }

  nPoints_ = areas.size();
  free_Gauss_Square(&Q,g+1);  
}


void WaveletBEM::constructSystemMatrix(GreensFunctionInterface &iface){

  generate_elementlist(&E_,P_,F_,p_,M_);
  generate_waveletlist(&W_,E_,p_,M_);
  set_quadrature_level(W_,E_,p_,M_);
  simplify_waveletlist(W_,E_,p_,M_);
  complete_elementlist(W_,E_,p_,M_);
  
  // Default GFs ok for the internal part
  apriori1_ = compression(&S_i_,W_,E_,p_,M_);
  WEM(&S_i_,W_,E_,T_,p_,M_,SingleLayerInt,DoubleLayerInt,2*M_PI);
  aposteriori1_ = postproc(&S_i_,W_,E_,p_,M_);

  // Use interface for the external part
  gf = &iface;
  apriori2_ = compression(&S_e_,W_,E_,p_,M_);
  WEM(&S_e_,W_,E_,T_,p_,M_,SingleLayer,DoubleLayer,-2*M_PI);
  aposteriori2_ = postproc(&S_e_,W_,E_,p_,M_);
  
  systemMatricesInitialized_ = true;
}


void WaveletBEM::solveForPotential(VectorXd &potential, VectorXd &charges){

  double *rhs;
  double *u = (double*) calloc(nf_,sizeof(double));
  double *v = (double*) calloc(nf_,sizeof(double));

  WEMRHS2M(&rhs,W_,E_,T_,p_,M_,potential.data(),quadratureLevel_);
  int iters = WEMPCG(&S_i_,rhs,u,eps_,p_,M_);
  memset(rhs,0,nf_*sizeof(double));
  for(unsigned int i=0; i<nf_; i++)
    for(unsigned int j=0; j<S_e_.row_number[i]; j++) 
      rhs[i] += S_e_.value1[i][j] * u[S_e_.index[i][j]];

  iters = WEMPGMRES3(&S_i_,&S_e_,rhs,v,eps_,p_,M_);
  for(unsigned int i=0; i<nf_; i++) u[i] -= 4*M_PI*v[i];
  tdwtKon(u,M_,nf_);

  // Interpolate charges

  cubature *Q;
  init_Gauss_Square(&Q,quadratureLevel_+1);

  vector2 s;
  vector2 t;
  int index = 0;
  int zi = 0;
  double h = 1.0/(1<<M_);

  for (unsigned int i1=0; i1<p_; i1++)
    {  s.y = 0;
      for (int i2=0; i2<(1<<M_); i2++)
	{  s.x = 0;
	  for (int i3=0; i3<(1<<M_); i3++)
	    {  for (unsigned int k=0; k<Q[quadratureLevel_].nop; k++)
		{  
		  t = vector2_add(s,vector2_Smul(h,Q[quadratureLevel_].xi[k]));
		  charges(index) = Q[quadratureLevel_].w[k]*u[zi] * h;
		  index ++;
		}
	      s.x += h;
	      zi++;
	    }
	  s.y += h;
	}
    }
    
  free_Gauss_Square(&Q,quadratureLevel_+1);
  free(rhs);
  free(u);
  free(v);

}

void WaveletBEM::writeChargesToDisk(VectorXd &charges){
  if((uint)charges.rows() != nf_*(quadratureLevel_+1)*(quadratureLevel_+1)){
    std::cerr << "WaveletBEM::writeChargesToDisk: Charges do" 
	      << " not match the operator definition." << std::endl;
    std::cerr << charges.rows() << " " << nf_ << " " << quadratureLevel_+1 << "\n";
    return;
  }

  std::ofstream file;
  file.open("charges.out",std::ios_base::out);

  file << nf_ << std::endl;

  uint p = 0;
  for(uint i = 0; i < nf_; i++){
    double avg = 0.0;
    for(uint j = 0; j < (quadratureLevel_+1)*(quadratureLevel_+1); j++){
      avg += charges(p);
      p++;
    }
    avg /= (quadratureLevel_+1)*(quadratureLevel_+1);
    file << avg << std::endl;
  }

  file.close();
}



void WaveletBEM::printInfo(std::ostream &out){
  out << " Operator parameters:" << std::endl
      << "   Patch level      : " << M_ << std::endl
      << "   Quadrature points: " << nPoints_ << std::endl
      << std::endl
      << "   Parameter a      : " << a << std::endl
      << "   Parameter b      : " << b << std::endl
      << "   Parameter dp     : " << dp << std::endl
      << std::endl
      << " System matrix sparsities:" << std::endl
      << "   A priori, matrix 1     : " << apriori1_ << std::endl
      << "   A posteriori, matrix 1 : " << aposteriori1_ << std::endl
      << "   A priori, matrix 2     : " << apriori2_ << std::endl
      << "   A posteriori, matrix 2 : " << aposteriori2_ << std::endl;
}

