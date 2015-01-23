
#include "LinAnsatzFunction.hpp"
#include "Compression.hpp"
#include "Transformations.hpp"
#include "Mask.hpp"
#include "Vector3.hpp"
#include "Vector2.hpp"

#include <ostream>
#include <cstdio>
#include <cstdlib>
#include "string.h"

// Ansatzfunktion 0
double Phi0(Vector2 a) {
	return ((1 - a.x) * (1 - a.y));
}
double dPhi0_dx(Vector2 a) {
	return (a.y - 1);
}
double dPhi0_dy(Vector2 a) {
	return (a.x - 1);
}

// Ansatzfunktion 1
double Phi1(Vector2 a) {
	return (a.x * (1 - a.y));
}
double dPhi1_dx(Vector2 a) {
	return (1 - a.y);
}
double dPhi1_dy(Vector2 a) {
	return (-a.x);
}

// Ansatzfunktion 2
double Phi2(Vector2 a) {
	return (a.x * a.y);
}
double dPhi2_dx(Vector2 a) {
	return (a.y);
}
double dPhi2_dy(Vector2 a) {
	return (a.x);
}

// Ansatzfunktion 3
double Phi3(Vector2 a) {
	return ((1 - a.x) * a.y);
}
double dPhi3_dx(Vector2 a) {
	return (-a.y);
}
double dPhi3_dy(Vector2 a) {
	return (1 - a.x);
}

// updated c_{i,j} by weight*phi_i(xi)*phi_j(eta)
void Phi_times_Phi(double *c, double weight, Vector2 xi, Vector2 eta){
  double a0, a1, a2, a3;
	double b0, b1, b2, b3;

	a0 = weight * (1 - xi.x) * (1 - xi.y);
	a1 = weight * xi.x * (1 - xi.y);
	a2 = weight * xi.x * xi.y;
	a3 = weight * (1 - xi.x) * xi.y;

	b0 = (1 - eta.x) * (1 - eta.y);
	b1 = eta.x * (1 - eta.y);
	b2 = eta.x * eta.y;
	b3 = (1 - eta.x) * eta.y;

	c[0] += a0 * b0;
	c[1] += a0 * b1;
	c[2] += a0 * b2;
	c[3] += a0 * b3;
	c[4] += a1 * b0;
	c[5] += a1 * b1;
	c[6] += a1 * b2;
	c[7] += a1 * b3;
	c[8] += a2 * b0;
	c[9] += a2 * b1;
	c[10] += a2 * b2;
	c[11] += a2 * b3;
	c[12] += a3 * b0;
	c[13] += a3 * b1;
	c[14] += a3 * b2;
	c[15] += a3 * b3;
	return;
}

// updated c_{i,j} by weight * < curl[phi_i(xi)],curl[phi_j(eta)] >
void Curl_Phi_times_Curl_Phi(double *c, double weight, Vector2 xi, Vector2 eta,
		Vector3 dChi_dx_s, Vector3 dChi_dy_s, Vector3 dChi_dx_t,
		Vector3 dChi_dy_t)
{
	double a0, a1, a2, a3;
	double b0, b1, b2, b3;
	double c0, c1, c2, c3;

	a0 = weight * vector3Dot(dChi_dy_s, dChi_dy_t);
	a1 = weight * vector3Dot(dChi_dy_s, dChi_dx_t);
	a2 = weight * vector3Dot(dChi_dx_s, dChi_dy_t);
	a3 = weight * vector3Dot(dChi_dx_s, dChi_dx_t);

	b0 = (1 - eta.y) * a0 - (1 - eta.x) * a1;
	b1 = (1 - eta.y) * a0 + eta.x * a1;
	b2 = eta.y * a0 - eta.x * a1;
	b3 = eta.y * a0 + (1 - eta.x) * a1;

	c0 = (1 - eta.y) * a2 - (1 - eta.x) * a3;
	c1 = (1 - eta.y) * a2 + eta.x * a3;
	c2 = eta.y * a2 - eta.x * a3;
	c3 = eta.y * a2 + (1 - eta.x) * a3;

	c[0] += +(1 - xi.y) * b0 - (1 - xi.x) * c0;
	c[1] += -(1 - xi.y) * b1 + (1 - xi.x) * c1;
	c[2] += -(1 - xi.y) * b2 + (1 - xi.x) * c2;
	c[3] += +(1 - xi.y) * b3 - (1 - xi.x) * c3;

	c[4] += -(1 - xi.y) * b0 - xi.x * c0;
	c[5] += +(1 - xi.y) * b1 + xi.x * c1;
	c[6] += +(1 - xi.y) * b2 + xi.x * c2;
	c[7] += -(1 - xi.y) * b3 - xi.x * c3;

	c[8] += -xi.y * b0 + xi.x * c0;
	c[9] += +xi.y * b1 - xi.x * c1;
	c[10] += +xi.y * b2 - xi.x * c2;
	c[11] += -xi.y * b3 + xi.x * c3;

	c[12] += +xi.y * b0 + (1 - xi.x) * c0;
	c[13] += -xi.y * b1 - (1 - xi.x) * c1;
	c[14] += -xi.y * b2 - (1 - xi.x) * c2;
	c[15] += +xi.y * b3 + (1 - xi.x) * c3;
	return;
}

LinAnsatzFunction :: LinAnsatzFunction(){
	nLevels = 0;
	nFunctions = 0;
	nPatches = 0;
  /// @todo check that the grade in the Interpolation is correct
  // now using grade = 4!!! (coefficient is the log2 (grade)
	//interCoeff  = new Interpolation(pPointsIn, 2, NEWTON, nLevels, nPatches);
	noPhi = 4;

  interCoeff = NULL;

  B = NULL;

  dp = 2.25;
  td = 4;
  a = 1.25; ///< compression constant,  a > 1
  b = 0.001; ///< compression constant, 0 < b < 1

  quadratureLevel_=2;
  //quadratureLevel_=1;
  G = (SparseMatrix*) malloc(sizeof(SparseMatrix));
}

LinAnsatzFunction :: LinAnsatzFunction(const Compression & _comp){
	nLevels = 0;
	nFunctions = 0;
	nPatches = 0;
  /// @todo check that the grade in the Interpolation is correct
  // now using grade = 4!!! (coefficient is the log2 (grade)
	//interCoeff  = new Interpolation(pPointsIn, 2, NEWTON, nLevels, nPatches);
	noPhi = 4;

  interCoeff = NULL;

  B = NULL;

  td = 4;
  dp = _comp.aPrioridPrime;
  a = _comp.aPrioriA; ///< compression constant,  a > 1
  b = _comp.aPosterioriB; ///< compression constant, 0 < b < 1

  quadratureLevel_=2;
  //quadratureLevel_=1;
  G = (SparseMatrix*) malloc(sizeof(SparseMatrix));
}


/// constructor of linear basis functions
LinAnsatzFunction :: LinAnsatzFunction(unsigned int _p, unsigned int _m, unsigned int _nf, double _a, double _b, double _dp, Vector3 *** pPointsIn){
	nLevels = _m;
	nFunctions = _nf;
	nPatches = _p;
  /// @todo check that the grade in the Interpolation is correct
  // now using grade = 4!!! (coefficient is the log2 (grade)
	interCoeff  = new Interpolation(pPointsIn, 2, NEWTON, nLevels, nPatches);
	noPhi = 4;


  B = NULL;

  dp = _dp;
  td = 4;
  a = _a; ///< compression constant,  a > 1
  b = _b; ///< compression constant, 0 < b < 1

  //quadratureLevel_=1;
  quadratureLevel_=2;
  G = (SparseMatrix*) malloc(sizeof(SparseMatrix));
}

/**
 * helper function that stores into the specific vector of the AnsatzFunction
 * the values of the integral
 **/
void LinAnsatzFunction::calculateIntegral(double* a, double *c){
  unsigned int k;
  double e[4];
  double m;
  unsigned int sizeC = 3*noPhi*noPhi;
  for( k = 0; k < sizeC; k+=4){
    m = 0.25 * (a[k+2]+a[sizeC+k+3]+a[2*sizeC+k] + a[3*sizeC+k+1]);
    e[0] = 0.5*(a[        k+1]+a[  sizeC+k]);
    e[1] = 0.5*(a[  sizeC+k+2]+a[2*sizeC+k+1]);
    e[2] = 0.5*(a[2*sizeC+k+3]+a[3*sizeC+k+2]);
    e[3] = 0.5*(a[3*sizeC+k  ]+a[        k+3]);

    c[k  ] = 0.5*(a[        k  ]+e[3]+e[0]+m);
    c[k+1] = 0.5*(a[1*sizeC+k+1]+e[0]+e[1]+m);
    c[k+2] = 0.5*(a[2*sizeC+k+2]+e[1]+e[2]+m);
    c[k+3] = 0.5*(a[3*sizeC+k+3]+e[2]+e[3]+m);
  }
}

/**
 * helper function that adds to the vector c the weight times the basis
 * functions evaluated at point xi
 */
void LinAnsatzFunction::calculateCRHS(double* c, double w, Vector2 xi){
  c[0] += w * Phi0(xi);
  c[1] += w * Phi1(xi);
  c[2] += w * Phi2(xi);
  c[3] += w * Phi3(xi);
  return;
}

/**
 * helper function that stores in the matrix y the influences of the children
 * elements on the value at the point i
 */
void LinAnsatzFunction::calculateYRHS(double** y, int i){
  double e[4];
  double m;
    
  m = 0.25 * (y[elementTree.element[i].son[0]][2]+y[elementTree.element[i].son[1]][3]+y[elementTree.element[i].son[2]][0]+y[elementTree.element[i].son[3]][1]);
  e[0] = 0.5*(y[elementTree.element[i].son[0]][1] + y[elementTree.element[i].son[1]][0]);
  e[1] = 0.5*(y[elementTree.element[i].son[1]][2] + y[elementTree.element[i].son[2]][1]);
  e[2] = 0.5*(y[elementTree.element[i].son[2]][3] + y[elementTree.element[i].son[3]][2]);
  e[3] = 0.5*(y[elementTree.element[i].son[3]][0] + y[elementTree.element[i].son[0]][3]);
    
  y[i][0] = 0.5*(y[elementTree.element[i].son[0]][0]+e[3]+e[0]+m);
  y[i][1] = 0.5*(y[elementTree.element[i].son[1]][1]+e[0]+e[1]+m);
  y[i][2] = 0.5*(y[elementTree.element[i].son[2]][2]+e[1]+e[2]+m);
  y[i][3] = 0.5*(y[elementTree.element[i].son[3]][3]+e[2]+e[3]+m);

  //printf("%lf %lf %lf %lf %d %d\n",y[i][0], y[i][1], y[i][2], y[i][3],elementTree.element[i].son[0], i);
}

/**
 * helper function that calculates the energy using the values in u and the
 * basis function evaluations in point xi
 */
double LinAnsatzFunction::calculateUEnergy(double *u, Vector2 xi, unsigned int zi){
  double U = 0.0;

  const unsigned int displacement = nPatches*((1<<(2*nLevels))-1)/3;
  et_node *pF = &elementTree.element[displacement+zi];
  U = u[pF->vertex[0]] *Phi0(xi)
    + u[pF->vertex[1]] *Phi1(xi)
    + u[pF->vertex[2]] *Phi2(xi)
    + u[pF->vertex[3]] *Phi3(xi);
  return U;
}


/**
 * function that creates the gramian matrix <phi_i, phi_j>
 */
void LinAnsatzFunction::createGram(unsigned int size, unsigned int maxRowNum){
  initSparse(G,size,size,maxRowNum);
  
  //do singleScaleGramian
  unsigned int n = 1<<nLevels;
  const unsigned int displacement = nPatches*(n*n-1)/3;
  double c[16];
  
  // initialization
  c[0] = c[5] = c[10] = c[15] = 1./9;
  c[1] = c[4] = c[ 9] = c[12] = 1./18;
  c[2] = c[7] = c[ 8] = c[13] = 1./36;
  c[3] = c[6] = c[11] = c[14] = 1./18;

  // construction of mass matrix
  for (unsigned int i=0; i<nPatches*n*n; ++i){
    for (unsigned int j=0; j<16; ++j)
      // L^2 normalized entries
      addSparse(G,elementTree.element[displacement+i].vertex[j%4],elementTree.element[displacement+i].vertex[j/4],c[j]);
  }
  return;
}

// releases the memory of the gramian matrix
void LinAnsatzFunction::freeGram(){
  freeSparse(G);
}

/// inverse wavelet transform
void LinAnsatzFunction :: tdwt(double *a){
  unsigned int ***C;
  unsigned int n, S, arg;

  SparseMatrix T, L;
  double *b;

  // loop over grid
  generateTopology(&C);
  //multiple(&C,&Z);
  b = (double*) malloc(waveletList.sizeWaveletList*sizeof(double));
  for(unsigned int m = minLevel; m <= nLevels; ++m){
    dwtMask(&T, &L, m, nLevels, this);
    n = 1<<m;
    S = 1<<(nLevels-m);
    
    // for scaling function set b = 0
    for(unsigned int i1 = 0; i1< nPatches; ++i1)
      for(unsigned int i2 = 0; i2 <= n; ++i2)
        for(unsigned int i3 = 0; i3 <= n; ++i3)
          b[C[i1][S*i2][S*i3]] = 0;

    // build scalarfunction and mother wavelets
    for(unsigned int i1 = 0; i1< nPatches; ++i1){
      for(unsigned int i2 = 0; i2 <= n; ++i2){
        for(unsigned int i3 = 0; i3 <= n; ++i3){
          // build tensor products psi_m(s) * phi_m+1 (t)
          if(i3%2 ==1){
            for(unsigned int s = 0; s < T.row_number[i3]; ++s){
              arg = C[i1][S*i2][S*T.index[i3][s]];
              b[arg] += T.value[i3][s] * a[C[i1][S*i2][S*i3]];
            }
          } else {
            // build tensor products psi_m(s) * phi_m (t) and phi_m (s) * phi_m (t)
            for(unsigned int t = 0; t < T.row_number[i2];++t){
              for(unsigned int s = 0; s < T.row_number[i3]; ++s){
                arg = C[i1][S*T.index[i2][t]][S*T.index[i3][s]];
                b[arg] += T.value[i3][s]*T.value[i2][t]*a[C[i1][S*i2][S*i3]];
              }
            }
          }
        }
      }
    }
    // update a
    for(unsigned int i1 = 0; i1< nPatches; ++i1)
      for(unsigned int i2 = 0; i2 <= n; ++i2)
        for(unsigned int i3 = 0; i3 <= n; ++i3)
          a[C[i1][S*i2][S*i3]] = 0.5*b[C[i1][S*i2][S*i3]];

    freeSparse(&T);
  }
  free(C);
  free(b);
  return;
}

/// wavelet transform
void LinAnsatzFunction :: dwt(double *a){
  unsigned int ***C;
  unsigned int n, S, arg;

  SparseMatrix T, L;
  double *b;

  // loop over grid
  //multiple(&C,&Z);
  generateTopology(&C);
  b = (double*) calloc(waveletList.sizeWaveletList,sizeof(double));
  for(unsigned int m = nLevels; m >= minLevel; --m){
    dwtMask(&T, &L, m, nLevels, this);
    n = 1<<m;
    S = 1<<(nLevels-m);

    // build scalarfunction and mother wavelets
    for(unsigned int i1 = 0; i1< nPatches; ++i1){
      for(unsigned int i2 = 0; i2 <= n; ++i2){
        for(unsigned int i3 = 0; i3 <= n; ++i3){
          // build tensor products psi_m(s) * phi_m+1 (t)
          if(i3%2 ==1){
            for(unsigned int s = 0; s < T.row_number[i3]; ++s){
              arg = C[i1][S*i2][S*T.index[i3][s]];
              b[C[i1][S*i2][S*i3]] += T.value[i3][s] * a[arg];
                //printf("DWT1: %d %lf %d %lf %lf\n", arg,T.value[i3][s],C[i1][S*i2][S*i3], b[C[i1][S*i2][S*i3]], a[arg]);
            }
          } else {
            // build tensor products psi_m(s) * phi_m (t) and phi_m (s) * phi_m (t)
            for(unsigned int t = 0; t < T.row_number[i2];++t){
              for(unsigned int s = 0; s < T.row_number[i3]; ++s){
                arg = C[i1][S*T.index[i2][t]][S*T.index[i3][s]];
                //printf("DWT: %d %lf %lf %d %lf %lf\n", arg,T.value[i3][s],T.value[i2][t],C[i1][S*i2][S*i3], b[C[i1][S*i2][S*i3]], a[arg]);
                b[C[i1][S*i2][S*i3]] += T.value[i3][s]*T.value[i2][t]*a[arg];
                //printf("DWT: %d %lf %lf %d %lf %lf\n", arg,T.value[i3][s],T.value[i2][t],C[i1][S*i2][S*i3], b[C[i1][S*i2][S*i3]], a[arg]);
              }
            }
          }
        }
      }
    }
    // update a
    for(unsigned int i1 = 0; i1< nPatches; ++i1)
      for(unsigned int i2 = 0; i2 <= n; ++i2)
        for(unsigned int i3 = 0; i3 <= n; ++i3)
          a[C[i1][S*i2][S*i3]] = 0.5*b[C[i1][S*i2][S*i3]];

    // for scaling function set b = 0
    for(unsigned int i1 = 0; i1< nPatches; ++i1)
      for(unsigned int i2 = 0; i2 <= n; i2+=2)
        for(unsigned int i3 = 0; i3 <= n; i3+=2)
          b[C[i1][S*i2][S*i3]] = 0;
    freeSparse(&T);
  }
  free(C);
  free(b);
  return;
}

/**
 * helper function permutating the integral according to the AnsatzFunction used
 **/
void LinAnsatzFunction::permutate(double *a, double *b){
  a[ 0] = b[ 0];
  a[ 1] = b[ 4];
  a[ 2] = b[ 8];
  a[ 3] = b[12];
  a[ 4] = b[ 1];
  a[ 5] = b[ 5];
  a[ 6] = b[ 9];
  a[ 7] = b[13];
  a[ 8] = b[ 2];
  a[ 9] = b[ 6];
  a[10] = b[10];
  a[11] = b[14];
  a[12] = b[ 3];
  a[13] = b[ 7];
  a[14] = b[11];
  a[15] = b[15];

  a[16] = b[32];
  a[17] = b[36];
  a[18] = b[40];
  a[19] = b[44];
  a[20] = b[33];
  a[21] = b[37];
  a[22] = b[41];
  a[23] = b[45];
  a[24] = b[34];
  a[25] = b[38];
  a[26] = b[42];
  a[27] = b[46];
  a[28] = b[35];
  a[29] = b[39];
  a[30] = b[43];
  a[31] = b[47];

  a[32] = b[16];
  a[33] = b[20];
  a[34] = b[24];
  a[35] = b[28];
  a[36] = b[17];
  a[37] = b[21];
  a[38] = b[25];
  a[39] = b[29];
  a[40] = b[18];
  a[41] = b[22];
  a[42] = b[26];
  a[43] = b[30];
  a[44] = b[19];
  a[45] = b[23];
  a[46] = b[27];
  a[47] = b[31];
}

/**
 * calculates the quadrature grade for elements with distance dist and levels
 * level1 and level2 using precition 2^(-alpha*level)
 *
 * @param dist distance between elements
 * @param alpha precision parameter
 * @param level1 level element1
 * @param level2 level element2
 * @param g1 required quadrature grade
 * @param g2 required quadrature grade
 */
void LinAnsatzFunction::quadratureGrade(signed int *g1, signed int*g2, int level1, int level2, double dist, double alpha) {
	dist = (dist < 1) ? log(dist) / log(2) : 0;

	// quadrature grade g1
	*g1 = (signed int) (-0.5*(alpha + level2)/(dist+level1+2));
	if(*g1 < 1) *g1 = 1;

	// quadrature grade g1
	*g2 = (signed int) (-0.5*(alpha + level1)/(dist+level2+2));
	if(*g2 < 1) *g2 = 1;

	return;
}

/// integration function for - no problem quadrature - modified scalar product
void LinAnsatzFunction::integrateNoProblem(double *c, unsigned int i1, unsigned int i2, unsigned int g1, Cubature *Q1, Cubature *Q2, double SingleLayer(Vector3, Vector3), double DoubleLayer(Vector3, Vector3, Vector3)) {
	double d1, d2, d3, hs;
	Vector2 xi, eta;
	Vector3 y, n_y;
  double ht = 1./(1<<elementTree.element[i2].level);
  unsigned int sizeC = 3*noPhi*noPhi;

  //calculate border integrals
  if(pRandWerte[g1].noP == 0){
    pRandWerte[g1].noP = Q1->noP;
    hs = 1./(1<<elementTree.element[i1].level);
    for(unsigned int i = 0; i < Q1->noP; ++i){
      xi.x = hs*(elementTree.element[i1].index_s+Q1->xi[i].x);
      xi.y = hs*(elementTree.element[i1].index_t+Q1->xi[i].y);

      pRandWerte[g1].Chi[i] = interCoeff->Chi(xi, elementTree.element[i1].patch);
      pRandWerte[g1].n_Chi[i] = interCoeff->n_Chi(xi, elementTree.element[i1].patch);
      pRandWerte[g1].det_dChi[i] = hs * Q1->weight[i];
    }
  }

	memset(c, 0, sizeC * sizeof(double));
	for (unsigned int i = 0; i < Q2->noP; ++i) {
    eta.x = ht * (elementTree.element[i2].index_s + Q2->xi[i].x);
    eta.y = ht * (elementTree.element[i2].index_t + Q2->xi[i].y);
    y = interCoeff->Chi(eta, elementTree.element[i2].patch);
    n_y = interCoeff->n_Chi(eta, elementTree.element[i2].patch);

		for (unsigned int j = 0; j < pRandWerte[g1].noP; ++j) {
			d1 = ht * Q2->weight[i] *pRandWerte[g1].det_dChi[j] * SingleLayer(pRandWerte[g1].Chi[j], y);
			d2 = ht * Q2->weight[i] *pRandWerte[g1].det_dChi[j] * DoubleLayer(pRandWerte[g1].Chi[j], y, n_y);
			d3 = ht * Q2->weight[i] *pRandWerte[g1].det_dChi[j] * DoubleLayer(y, pRandWerte[g1].Chi[j], pRandWerte[g1].n_Chi[j]);
			Phi_times_Phi(&c[ 0], d1, Q1->xi[j], Q2->xi[i]);
			Phi_times_Phi(&c[16], d2, Q1->xi[j], Q2->xi[i]);
			Phi_times_Phi(&c[32], d3, Q1->xi[j], Q2->xi[i]);
		}
	}
	return;
}

/// integration function for - same patches
void LinAnsatzFunction::integratePatch(double *c, unsigned int i1, Cubature * Q, double SingleLayer(Vector3, Vector3), double DoubleLayer(Vector3, Vector3, Vector3), double Identity) {
	double d1, d2, d3, w;
	double t1, t2, t3, t4;
	Vector2 s, xi, eta, a, b;
	Vector3 x, y;
	double h = 1. / (1 << elementTree.element[i1].level);
  unsigned int sizeC = 3*noPhi*noPhi;

	memset(c, 0, sizeC * sizeof(double));
  s = Vector2(h*elementTree.element[i1].index_s, h*elementTree.element[i1].index_t);

	for (unsigned int i = 0; i < Q->noP; ++i) {
		xi = Q->xi[i];
		w = h * h * Q->weight[i] * xi.x * (1 - xi.x) * (1 - xi.x * xi.y);
		for (unsigned int j = 0; j < Q->noP; ++j) {
			eta = Q->xi[j];
			t1 = eta.x * (1 - xi.x);
			t2 = eta.y * (1 - xi.x * xi.y);
			t3 = t1 + xi.x;
			t4 = t2 + xi.x * xi.y;

			a.x = s.x + h * t1;
			a.y = s.y + h * t2;
			b.x = s.x + h * t3;
			b.y = s.y + h * t4;
			x = interCoeff->Chi(a, elementTree.element[i1].patch);
			y = interCoeff->Chi(b, elementTree.element[i1].patch);
			d1 = w * Q->weight[j] * SingleLayer(x, y);
			d2 = w * Q->weight[j] * DoubleLayer(x, y, interCoeff->n_Chi(b, elementTree.element[i1].patch));
			d3 = w * Q->weight[j] * DoubleLayer(y, x, interCoeff->n_Chi(a, elementTree.element[i1].patch));
			Phi_times_Phi(&c[0], d1, Vector2(t1, t2), Vector2(t3, t4));
			Phi_times_Phi(&c[0], d1, Vector2(t3, t4), Vector2(t1, t2));
			Phi_times_Phi(&c[16], d2, Vector2(t1, t2), Vector2(t3, t4));
			Phi_times_Phi(&c[16], d3, Vector2(t3, t4), Vector2(t1, t2));

			a.y = s.y + h * t4;
			b.y = s.y + h * t2;
			x = interCoeff->Chi(a, elementTree.element[i1].patch);
			y = interCoeff->Chi(b, elementTree.element[i1].patch);
			d1 = w * Q->weight[j] * SingleLayer(x, y);
			d2 = w * Q->weight[j] * DoubleLayer(x, y, interCoeff->n_Chi(b, elementTree.element[i1].patch));
			d3 = w * Q->weight[j] * DoubleLayer(y, x, interCoeff->n_Chi(a, elementTree.element[i1].patch));
			Phi_times_Phi(&c[0], d1, Vector2(t1, t4), Vector2(t3, t2));
			Phi_times_Phi(&c[0], d1, Vector2(t3, t2), Vector2(t1, t4));
			Phi_times_Phi(&c[16], d2, Vector2(t1, t4), Vector2(t3, t2));
			Phi_times_Phi(&c[16], d3, Vector2(t3, t2), Vector2(t1, t4));

			a.x = s.x + h * t2;
			a.y = s.y + h * t1;
			b.x = s.x + h * t4;
			b.y = s.y + h * t3;
			x = interCoeff->Chi(a, elementTree.element[i1].patch);
			y = interCoeff->Chi(b, elementTree.element[i1].patch);
			d1 = w * Q->weight[j] * SingleLayer(x, y);
			d2 = w * Q->weight[j] * DoubleLayer(x, y, interCoeff->n_Chi(b, elementTree.element[i1].patch));
			d3 = w * Q->weight[j] * DoubleLayer(y, x, interCoeff->n_Chi(a, elementTree.element[i1].patch));
			Phi_times_Phi(&c[0], d1, Vector2(t2, t1), Vector2(t4, t3));
			Phi_times_Phi(&c[0], d1, Vector2(t4, t3), Vector2(t2, t1));
			Phi_times_Phi(&c[16], d2, Vector2(t2, t1), Vector2(t4, t3));
			Phi_times_Phi(&c[16], d3, Vector2(t4, t3), Vector2(t2, t1));

			a.y = s.y + h * t3;
			b.y = s.y + h * t1;
			x = interCoeff->Chi(a, elementTree.element[i1].patch);
			y = interCoeff->Chi(b, elementTree.element[i1].patch);
			d1 = w * Q->weight[j] * SingleLayer(x, y);
			d2 = w * Q->weight[j] * DoubleLayer(x, y, interCoeff->n_Chi(b, elementTree.element[i1].patch));
			d3 = w * Q->weight[j] * DoubleLayer(y, x, interCoeff->n_Chi(a, elementTree.element[i1].patch));
			Phi_times_Phi(&c[0], d1, Vector2(t2, t3), Vector2(t4, t1));
			Phi_times_Phi(&c[0], d1, Vector2(t4, t1), Vector2(t2, t3));
			Phi_times_Phi(&c[16], d2, Vector2(t2, t3), Vector2(t4, t1));
			Phi_times_Phi(&c[16], d3, Vector2(t4, t1), Vector2(t2, t3));
		}
	}

  w = Identity;
	c[16] += w/9;
	c[17] += w/18;
	c[18] += w/36;
	c[19] += w/18;
	c[20] += w/18;
	c[21] += w/9;
	c[22] += w/18;
	c[23] += w/36;
	c[24] += w/36;
	c[25] += w/18;
	c[26] += w/9;
	c[27] += w/18;
	c[28] += w/18;
	c[29] += w/36;
	c[30] += w/18;
	c[31] += w/9;

	// transposed entry
	c[32] = c[16];
	c[33] = c[20];
	c[34] = c[24];
	c[35] = c[28];
	c[36] = c[17];
	c[37] = c[21];
	c[38] = c[25];
	c[39] = c[29];
	c[40] = c[18];
	c[41] = c[22];
	c[42] = c[26];
	c[43] = c[30];
	c[44] = c[19];
	c[45] = c[23];
	c[46] = c[27];
	c[47] = c[31];
	return;
}

/// integration function for - common edge
void LinAnsatzFunction::integrateEdge(double *c, unsigned int i1, unsigned int i2, unsigned int ind_s, unsigned int ind_t, Cubature *Q, double SingleLayer(Vector3, Vector3), double DoubleLayer(Vector3, Vector3, Vector3), double Identity) {
	double d1, d2, d3, w, t1, t2, t3, t4;
	Vector2 xi, eta, a, b, u, v, s, t;
	Vector3 x, y;
	double h = 1. / (1 << elementTree.element[i1].level);
  unsigned int sizeC=3*noPhi*noPhi;

	memset(c, 0, sizeC * sizeof(double));
  s = Vector2(h*elementTree.element[i1].index_s, h* elementTree.element[i1].index_t); 
  t = Vector2(h*elementTree.element[i2].index_s, h* elementTree.element[i2].index_t);

	for (unsigned int i = 0; i < Q->noP; ++i) {
		xi = Q->xi[i];
		w = h * h * xi.y * xi.y * Q->weight[i];
		t1 = xi.x * (1 - xi.y);
		t2 = (1 - xi.x) * (1 - xi.y);

		for (unsigned int j = 0; j < Q->noP; ++j) {
			eta = vector2SMul(xi.y, Q->xi[j]);
			t3 = xi.x * (1 - eta.x);
			t4 = (1 - xi.x) * (1 - eta.x);

			a = tau(t1, eta.x, ind_s);
			b = tau(t2, eta.y, ind_t);
			u = kappa(s, a, h);
			v = kappa(t, b, h);
			x = interCoeff->Chi(u, elementTree.element[i1].patch);
			y = interCoeff->Chi(v, elementTree.element[i2].patch);
			d1 = w * Q->weight[j] * (1 - xi.y) * SingleLayer(x, y);
			d2 = w * Q->weight[j] * (1 - xi.y)
					* DoubleLayer(x, y, interCoeff->n_Chi(v, elementTree.element[i2].patch));
			d3 = w * Q->weight[j] * (1 - xi.y)
					* DoubleLayer(y, x, interCoeff->n_Chi(u, elementTree.element[i1].patch));
			Phi_times_Phi(&c[ 0], d1, a, b);
			Phi_times_Phi(&c[16], d2, a, b);
			Phi_times_Phi(&c[32], d3, a, b);

			a = tau(1 - t1, eta.x, ind_s);
			b = tau(1 - t2, eta.y, ind_t);
			u = kappa(s, a, h);
			v = kappa(t, b, h);
			x = interCoeff->Chi(u, elementTree.element[i1].patch);
			y = interCoeff->Chi(v, elementTree.element[i2].patch);
			d1 = w * Q->weight[j] * (1 - xi.y) * SingleLayer(x, y);
			d2 = w * Q->weight[j] * (1 - xi.y)
					* DoubleLayer(x, y, interCoeff->n_Chi(v, elementTree.element[i2].patch));
			d3 = w * Q->weight[j] * (1 - xi.y)
					* DoubleLayer(y, x, interCoeff->n_Chi(u, elementTree.element[i1].patch));
			Phi_times_Phi(&c[0], d1, a, b);
			Phi_times_Phi(&c[16], d2, a, b);
			Phi_times_Phi(&c[32], d3, a, b);

			a = tau(t3, xi.y, ind_s);
			b = tau(t4, eta.y, ind_t);
			u = kappa(s, a, h);
			v = kappa(t, b, h);
			x = interCoeff->Chi(u, elementTree.element[i1].patch);
			y = interCoeff->Chi(v, elementTree.element[i2].patch);
			d1 = w * Q->weight[j] * (1 - eta.x) * SingleLayer(x, y);
			d2 = w * Q->weight[j] * (1 - eta.x)
					* DoubleLayer(x, y, interCoeff->n_Chi(v, elementTree.element[i2].patch));
			d3 = w * Q->weight[j] * (1 - eta.x)
					* DoubleLayer(y, x, interCoeff->n_Chi(u, elementTree.element[i1].patch));
			Phi_times_Phi(&c[ 0], d1, a, b);
			Phi_times_Phi(&c[16], d2, a, b);
			Phi_times_Phi(&c[32], d3, a, b);

			a = tau(1 - t3, xi.y, ind_s);
			b = tau(1 - t4, eta.y, ind_t);
			u = kappa(s, a, h);
			v = kappa(t, b, h);
			x = interCoeff->Chi(u, elementTree.element[i1].patch);
			y = interCoeff->Chi(v, elementTree.element[i2].patch);
			d1 = w * Q->weight[j] * (1 - eta.x) * SingleLayer(x, y);
			d2 = w * Q->weight[j] * (1 - eta.x)
					* DoubleLayer(x, y, interCoeff->n_Chi(v, elementTree.element[i2].patch));
			d3 = w * Q->weight[j] * (1 - eta.x)
					* DoubleLayer(y, x, interCoeff->n_Chi(u, elementTree.element[i1].patch));
			Phi_times_Phi(&c[ 0], d1, a, b);
			Phi_times_Phi(&c[16], d2, a, b);
			Phi_times_Phi(&c[32], d3, a, b);

			a = tau(t4, eta.y, ind_s);
			b = tau(t3, xi.y, ind_t);
			u = kappa(s, a, h);
			v = kappa(t, b, h);
			x = interCoeff->Chi(u, elementTree.element[i1].patch);
			y = interCoeff->Chi(v, elementTree.element[i2].patch);
			d1 = w * Q->weight[j] * (1 - eta.x) * SingleLayer(x, y);
			d2 = w * Q->weight[j] * (1 - eta.x)
					* DoubleLayer(x, y, interCoeff->n_Chi(v, elementTree.element[i2].patch));
			d3 = w * Q->weight[j] * (1 - eta.x)
					* DoubleLayer(y, x, interCoeff->n_Chi(u, elementTree.element[i1].patch));
			Phi_times_Phi(&c[ 0], d1, a, b);
			Phi_times_Phi(&c[16], d2, a, b);
			Phi_times_Phi(&c[32], d3, a, b);

			a = tau(1 - t4, eta.y, ind_s);
			b = tau(1 - t3, xi.y, ind_t);
			u = kappa(s, a, h);
			v = kappa(t, b, h);
			x = interCoeff->Chi(u, elementTree.element[i1].patch);
			y = interCoeff->Chi(v, elementTree.element[i2].patch);
			d1 = w * Q->weight[j] * (1 - eta.x) * SingleLayer(x, y);
			d2 = w * Q->weight[j] * (1 - eta.x)
					* DoubleLayer(x, y, interCoeff->n_Chi(v, elementTree.element[i2].patch));
			d3 = w * Q->weight[j] * (1 - eta.x)
					* DoubleLayer(y, x, interCoeff->n_Chi(u, elementTree.element[i1].patch));
			Phi_times_Phi(&c[ 0], d1, a, b);
			Phi_times_Phi(&c[16], d2, a, b);
			Phi_times_Phi(&c[32], d3, a, b);
		}
	}
	return;
}

/// integration function for - common node in origin
void LinAnsatzFunction::integratePoint(double *c, unsigned int i1, unsigned int i2, unsigned int ind_s,
		unsigned int ind_t, Cubature *Q, double SingleLayer(Vector3, Vector3), double DoubleLayer(Vector3, Vector3, Vector3), double Identity) {
	double d1, d2, d3, w;
	Vector2 xi, eta, a, u, a1, a2, b1, b2, s, t;
	Vector3 x1, n_x1, x2, n_x2, y1, n_y1, y2, n_y2, z, n_z;
	double h = 1. / (1 << elementTree.element[i1].level);
  unsigned int sizeC = 3*noPhi*noPhi;

	memset(c, 0, sizeC * sizeof(double));
  s = Vector2(h*elementTree.element[i1].index_s, h*elementTree.element[i1].index_t);
  t = Vector2(h*elementTree.element[i2].index_s, h*elementTree.element[i2].index_t);

	for (unsigned int i = 0; i < Q->noP; ++i) {
		xi = Q->xi[i];
		w = h * h * pow(xi.x, 3) * Q->weight[i];
		xi.y *= xi.x;
		a1 = tau(xi.x, xi.y, ind_s);
		a2 = tau(xi.y, xi.x, ind_s);
		b1 = tau(xi.x, xi.y, ind_t);
		b2 = tau(xi.y, xi.x, ind_t);

		u = kappa(s, a1, h);
		x1 = interCoeff->Chi(u, elementTree.element[i1].patch);
		n_x1 = interCoeff->n_Chi(u, elementTree.element[i1].patch);

		u = kappa(s, a2, h);
		x2 = interCoeff->Chi(u, elementTree.element[i1].patch);
		n_x2 = interCoeff->n_Chi(u, elementTree.element[i1].patch);

		u = kappa(t, b1, h);
		y1 = interCoeff->Chi(u, elementTree.element[i2].patch);
		n_y1 = interCoeff->n_Chi(u, elementTree.element[i2].patch);

		u = kappa(t, b2, h);
		y2 = interCoeff->Chi(u, elementTree.element[i2].patch);
		n_y2 = interCoeff->n_Chi(u, elementTree.element[i2].patch);

		for (unsigned int j = 0; j < Q->noP; ++j) {
			eta = vector2SMul(xi.x, Q->xi[j]);

			a = tau(eta.x, eta.y, ind_t);
			u = kappa(t, a, h);
			z = interCoeff->Chi(u, elementTree.element[i2].patch);
			n_z = interCoeff->n_Chi(u, elementTree.element[i2].patch);

			d1 = w * Q->weight[j] * SingleLayer(x1, z);
			d2 = w * Q->weight[j] * DoubleLayer(x1, z, n_z);
			d3 = w * Q->weight[j] * DoubleLayer(z, x1, n_x1);
			Phi_times_Phi(&c[ 0], d1, a1, a);
			Phi_times_Phi(&c[16], d2, a1, a);
			Phi_times_Phi(&c[32], d3, a1, a);

			d1 = w * Q->weight[j] * SingleLayer(x2, z);
			d2 = w * Q->weight[j] * DoubleLayer(x2, z, n_z);
			d3 = w * Q->weight[j] * DoubleLayer(z, x2, n_x2);
			Phi_times_Phi(&c[ 0], d1, a2, a);
			Phi_times_Phi(&c[16], d2, a2, a);
			Phi_times_Phi(&c[32], d3, a2, a);

			a = tau(eta.x, eta.y, ind_s);
			u = kappa(s, a, h);
			z = interCoeff->Chi(u, elementTree.element[i1].patch);
			n_z = interCoeff->n_Chi(u, elementTree.element[i1].patch);

			d1 = w * Q->weight[j] * SingleLayer(z, y1);
			d2 = w * Q->weight[j] * DoubleLayer(z, y1, n_y1);
			d3 = w * Q->weight[j] * DoubleLayer(y1, z, n_z);
			Phi_times_Phi(&c[ 0], d1, a, b1);
			Phi_times_Phi(&c[16], d2, a, b1);
			Phi_times_Phi(&c[32], d3, a, b1);

			d1 = w * Q->weight[j] * SingleLayer(z, y2);
			d2 = w * Q->weight[j] * DoubleLayer(z, y2, n_y2);
			d3 = w * Q->weight[j] * DoubleLayer(y2, z, n_z);
			Phi_times_Phi(&c[ 0], d1, a, b2);
			Phi_times_Phi(&c[16], d2, a, b2);
			Phi_times_Phi(&c[32], d3, a, b2);
		}
	}
	return;
}

void LinAnsatzFunction :: generateCanonicalSingleScaleBasis(Wavelet *W, unsigned int ***C, int m){
	unsigned int	n = 1 << m;                 // p*n*n elements on level m
	unsigned int	S = 1 << (nLevels-m); // step to next patch?
	unsigned int	ze;                         // element counter
	const unsigned int sizeNewWeights = 4;
	double *newWeights = (double*) calloc(sizeNewWeights,sizeof(double));

  // initialization
	for (unsigned int i1=0; i1<nPatches; ++i1) {
		for (unsigned int i2=0; i2<=n; ++i2) {
			for (unsigned int i3=0; i3<=n; ++i3) {
				W[C[i1][S*i2][S*i3]].noElements = 0;
				W[C[i1][S*i2][S*i3]].level = m;
			}
		}
	}
	
  // calculate singlescale basis on level m
	ze = nPatches*(n*n-1)/3; // zero element on level m
	for (unsigned int i1=0; i1<nPatches; ++i1) {
		for (unsigned int i2=0; i2<n; ++i2) {
			for (unsigned int i3=0; i3<n; ++i3) {
				//generate vector for weights
				newWeights[0] = 1;
				addElement(&W[C[i1][S* i2	][S* i3   ]],elementTree,newWeights,sizeNewWeights,ze);
				newWeights[0] = 0;
				newWeights[1] = 1;
				addElement(&W[C[i1][S* i2	][S*(i3+1)]],elementTree,newWeights,sizeNewWeights,ze);
				newWeights[1] = 0;
				newWeights[2] = 1;
				addElement(&W[C[i1][S*(i2+1)][S*(i3+1)]],elementTree,newWeights,sizeNewWeights,ze);
				newWeights[2] = 0;
				newWeights[3] = 1;
				addElement(&W[C[i1][S*(i2+1)][S* i3   ]],elementTree,newWeights,sizeNewWeights,ze);
				newWeights[3] = 0;
				++ze;
			}
		}
	}
  free(newWeights);
	return;
}

// refine coarse elements
void LinAnsatzFunction :: setQuadratureLevel() {
	unsigned int	ind;            // element under consideration
	unsigned int	minLevelLocal;  // minimal level for quadrature, composed of global value and number of refinements
	unsigned int	noe;            // number of elements of wavelet
	double			*w;               // weights of element functions

	minLevelLocal = minQuadratureLevel;// minimal level for quadrature, composed of global value and number of refinements
	if (minLevelLocal > nLevels) minLevelLocal = nLevels; // minimal level at most refinement level
	w = (double*) malloc(noPhi*sizeof(double));

	for (unsigned int i=0; (i<waveletList.sizeWaveletList)&&(waveletList.W[i].level<=minLevelLocal); ++i) {
		for (unsigned int j=waveletList.W[i].level; (j<=minLevelLocal); ++j) {
			noe = waveletList.W[i].noElements;
			for (unsigned int k=0; k<noe; ++k) {
        // element under consideration
				ind = waveletList.W[i].element[k];
        // is it too coarse?
				if (elementTree.element[ind].level < minLevelLocal) {
					waveletList.W[i].element[k] = elementTree.element[ind].son[0];
					memcpy(w,&waveletList.W[i].weight[k*noPhi],noPhi*sizeof(double));
					waveletList.W[i].element = (unsigned int*) realloc(waveletList.W[i].element,(waveletList.W[i].noElements+3)*sizeof(unsigned int));
					waveletList.W[i].weight  = (double*) realloc(waveletList.W[i].weight ,(waveletList.W[i].noElements+3)*noPhi*sizeof(double));

          // calculate new weights
					for (unsigned int l=1; l<4; ++l) {
						waveletList.W[i].element[waveletList.W[i].noElements] = elementTree.element[ind].son[l];
						waveletList.W[i].weight[waveletList.W[i].noElements*noPhi+ l-1	 ] = 0.25*(w[l-1]+w[l]);
						waveletList.W[i].weight[waveletList.W[i].noElements*noPhi+ l	 ] = 0.5*w[l]; 
						waveletList.W[i].weight[waveletList.W[i].noElements*noPhi+(l+1)%4] = 0.25*(w[l]+w[(l+1)%4]);
						waveletList.W[i].weight[waveletList.W[i].noElements*noPhi+(l+2)%4] = 0.125*(w[0]+w[1]+w[2]+w[3]);
						++waveletList.W[i].noElements;
					}
					waveletList.W[i].element[k] = elementTree.element[ind].son[0];
					waveletList.W[i].weight[k*noPhi+0] = 0.5*w[0];
					waveletList.W[i].weight[k*noPhi+1] = 0.25*(w[0]+w[1]);
					waveletList.W[i].weight[k*noPhi+2] = 0.125*(w[0]+w[1]+w[2]+w[3]);
					waveletList.W[i].weight[k*noPhi+3] = 0.25*(w[3]+w[0]);
				}
			}
		}
	}
#ifdef DEBUG
	FILE* debugFile = fopen("debug.out","a");
	fprintf(debugFile,">>> WAVELET_TREE_QUADRATURE\n");
	for(unsigned int m = 0; m<waveletList.sizeWaveletList; ++m){
		fprintf(debugFile,"%d %d %d %d\n", m, waveletList.W[m].level, waveletList.W[m].noElements, waveletList.W[m].noSons);
		for(unsigned int i1 = 0 ; i1< waveletList.W[m].noElements; ++i1){
			fprintf(debugFile,"%d %lf %lf %lf %lf ", waveletList.W[m].element[i1], waveletList.W[m].weight[i1*4], waveletList.W[m].weight[i1*4+1], waveletList.W[m].weight[i1*4+2], waveletList.W[m].weight[i1*4+3]);
		}
		for(unsigned int i1 = 0 ; i1< waveletList.W[m].noSons; ++i1){
			fprintf(debugFile,"%d ", waveletList.W[m].son[i1]);
		}
		fprintf(debugFile,"\n");
	}
	fprintf(debugFile,"<<< WAVELET_TREE_QUADRATURE\n");
	fclose(debugFile);
#endif

  // free auxiliary memory
  free(w);
	return;
}

void LinAnsatzFunction :: generateTopology(unsigned int ****C) {
	unsigned int n = 1<< nLevels;  // n*n elements on patch
	unsigned int ze = nPatches*(n*n-1)/3; // index in element list
#ifdef DEBUG
	FILE* debugFile = fopen("debug.out","a");
	fprintf(debugFile,">>> GENERATETOPOLOGY C\n");
#endif

  // allocate memory - done in one go
	(*C) = (unsigned int***) malloc(nPatches*sizeof(unsigned int*)+(nPatches*(n+1)*sizeof(unsigned int**))+(nPatches*(n+1)*(n+1)*sizeof(unsigned int)));
	for(unsigned int i1 = 0; i1 < nPatches; ++i1){
		(*C)[i1] = (unsigned int**)((*C)+nPatches)+i1*(n+1);
		for(unsigned int i2 = 0; i2 <= n; ++i2){
			(*C)[i1][i2] = (unsigned int*)((*C)+nPatches+nPatches*(n+1))+i1*(n+1)*(n+1)+i2*(n+1);
			memset((*C)[i1][i2],0,(n+1)*sizeof(unsigned int));
			for(unsigned int i3 = 0; i3 <= n; ++i3){
				if		((i2 <= n/2) && (i3 <= n/2)) (*C)[i1][i2][i3] = elementTree.element[ze].vertex[0];
				else if ((i2 <= n/2) && (i3 >  n/2)) (*C)[i1][i2][i3] = elementTree.element[ze].vertex[1];
				else if ((i2 >	n/2) && (i3 >  n/2)) (*C)[i1][i2][i3] = elementTree.element[ze].vertex[2];
				else								 (*C)[i1][i2][i3] = elementTree.element[ze].vertex[3];
				if (i3 != n/2) ++ze;
#ifdef DEBUG
				fprintf(debugFile, "%3d\t",(*C)[i1][i2][i3]);
#endif
			}
#ifdef DEBUG
			fprintf(debugFile, "\n");
#endif
			if (i2 == n/2) ze -= n;
		}
#ifdef DEBUG
		fprintf(debugFile, "\n");
#endif
	}
#ifdef DEBUG
	fprintf(debugFile,"<<< GENERATETOPOLOGY C\n");
	fclose(debugFile);
#endif
	return;
}

/// construction of WaveletList
void LinAnsatzFunction :: generateWaveletList() {
	unsigned int ***C;       // local basis list
	unsigned int n;          // n*n elements per patch on level m
	Wavelet *w;              // wavelet under consideration
	SparseMatrix Th, L;      // mask matrices
	unsigned int S;          // distance to next point
	unsigned int arg;        // temporary variable
	Wavelet *G;              // temporary wavelet list
	double *newWeights;      // new weights to be added to w
	unsigned int sizeNewWeights;// size of the weightsList
	
	sizeNewWeights = 4;
	newWeights = (double*) calloc(sizeNewWeights,sizeof(double));

	waveletList.sizeWaveletList = nNodes;
	// 1. initialization
	generateTopology(&C);
	G = (Wavelet*) calloc(waveletList.sizeWaveletList,sizeof(Wavelet));
	waveletList.W = (Wavelet*) calloc(waveletList.sizeWaveletList,sizeof(Wavelet));

	// 2. loop over grid levels
	for (int m=nLevels; m>=(int)minLevel; m--) {
    // calculate mask matrices Th and L
		dwtMask(&Th,&L,m,m,this);
		n = 1 << m;                 // nPatches*n*n elements on level m 
		S = 1 << (nLevels-m); // step size to next point

		// 2.1. build wavelets
		generateCanonicalSingleScaleBasis(G,C,m);
		for (unsigned int i1=0; i1<nPatches; ++i1) {
			for (unsigned int i2=0; i2<=n; ++i2) {
				for (unsigned int i3=0; i3<=n; ++i3){
          // wavelet under consideration
					w = &(waveletList.W)[C[i1][S*i2][S*i3]];
					w->level = m;			// level of wavelet 
          // add elements to wavelet
					if (i3%2 == 1) {
            // tensorproduct psi_{m}(s)*phi_{m+1}(t)
						for (unsigned int s=0; s<Th.row_number[i3]; ++s){
							arg = C[i1][S*i2][S*Th.index[i3][s]];
							//add_wavelet(w,&G[arg],E,0.5*Th.value[i3][s]);
							for(unsigned int l = 0; l < G[arg].noElements; ++l){
								memset(newWeights, 0, sizeNewWeights*sizeof(double));
								// generate the proper weights
								for(unsigned int k = 0; k < sizeNewWeights; ++k){
									newWeights[k] = 0.5*Th.value[i3][s]*G[arg].weight[l*sizeNewWeights+k];
								}
								addElement(w, elementTree, newWeights, sizeNewWeights, G[arg].element[l]);
							}
						}
					}
					else if (i2%2 == 1) {
            // tensorproducte phi_{m}(s)*psi_{m}(t)
						for (unsigned int t=0; t<Th.row_number[i2]; ++t) {
							for (unsigned int s=0; s<Th.row_number[i3]; ++s){
								arg = C[i1][S*Th.index[i2][t]][S*Th.index[i3][s]];
								//add_wavelet(w,&G[arg],E,0.5*Th.value[i3][s]*Th.value[i2][t]);
								for(unsigned int l = 0; l < G[arg].noElements; ++l){
									memset(newWeights, 0, sizeNewWeights*sizeof(double));
									// generate the proper weights
									for(unsigned int k = 0; k < sizeNewWeights; ++k){
										newWeights[k] = 0.5*Th.value[i3][s]*Th.value[i2][t]*G[arg].weight[l*sizeNewWeights+k];
									}
									addElement(w, elementTree, newWeights, sizeNewWeights, G[arg].element[l]);
								}
							}
						}
					}
          // calculate children
					if ((i3%4 == 2) && (i2%2 == 0)){
						if (i2 < n){
							w->noSons = 4;
							w->son = (unsigned int*) realloc(w->son,4*sizeof(unsigned int));
							w->son[0] = C[i1][S* i2   ][S*(i3-1)];
							w->son[1] = C[i1][S* i2   ][S*(i3+1)];
							w->son[2] = C[i1][S*(i2+1)][S*(i3-1)];
							w->son[3] = C[i1][S*(i2+1)][S*(i3+1)];
						}else{
							w->noSons = 2;
							w->son = (unsigned int*) realloc(w->son,2*sizeof(unsigned int));
							w->son[0] = C[i1][S*i2][S*(i3-1)];
							w->son[1] = C[i1][S*i2][S*(i3+1)];
						}
					} else if ((i2%4 == 2) && (i3%4 == 0)){
						if (i3 < n){
							w->noSons = 4;
							w->son = (unsigned int*) realloc(w->son,4*sizeof(unsigned int));
							w->son[0] = C[i1][S*(i2-1)][S* i3	];
							w->son[1] = C[i1][S*(i2-1)][S*(i3+2)];
							w->son[2] = C[i1][S*(i2+1)][S* i3	];
							w->son[3] = C[i1][S*(i2+1)][S*(i3+2)];
						} else {
							w->noSons = 2;
							w->son = (unsigned int*) realloc(w->son,2*sizeof(unsigned int));
							w->son[0] = C[i1][S*(i2-1)][S*i3];
							w->son[1] = C[i1][S*(i2+1)][S*i3];
						}
					}
				}
			}
		}
		freeSparse(&Th);
	}

  // 3. add scaling functions
	generateCanonicalSingleScaleBasis(waveletList.W,C,minLevel-1);

  // free memory
	for (unsigned int i1=0; i1<waveletList.sizeWaveletList; ++i1){
		free(G[i1].element);
		free(G[i1].weight);
		free(G[i1].son);
	}
  free(newWeights);
	free(G);
	free(C); // wird automatisch ge-freed am ende der funktion...
#ifdef DEBUG
	FILE* debugFile = fopen("debug.out","a");
	fprintf(debugFile,">>> WAVELET_TREE\n");
	for(unsigned int m = 0; m<waveletList.sizeWaveletList; ++m){
		fprintf(debugFile,"%d %d %d %d\n", m, waveletList.W[m].level, waveletList.W[m].noElements, waveletList.W[m].noSons);
		for(unsigned int i1 = 0 ; i1< waveletList.W[m].noElements; ++i1){
			fprintf(debugFile,"%d %lf %lf %lf %lf ", waveletList.W[m].element[i1], waveletList.W[m].weight[i1*4], waveletList.W[m].weight[i1*4+1], waveletList.W[m].weight[i1*4+2], waveletList.W[m].weight[i1*4+3]);
		}
		for(unsigned int i1 = 0 ; i1< waveletList.W[m].noSons; ++i1){
			fprintf(debugFile,"%d ", waveletList.W[m].son[i1]);
		}
		fprintf(debugFile,"\n");
	}
	fprintf(debugFile,"<<< WAVELET_TREE\n");
	fclose(debugFile);
#endif

	return;
}

/**
 * removes fine level wavelets that hold no new information - replaces them by
 * parent
 */
void LinAnsatzFunction :: simplifyWaveletList() {
	unsigned int	k;               // helper index
	signed int		j;               // index in wavelet/element list
	unsigned int	s1, s2, s3;      // indeces for children
	unsigned int	noe;             // number of elements in wavelet support
	unsigned int	*prototype;      // prototype index list
	unsigned int	prototype_number;// number of prototypes
	unsigned int	minLevelLocal;   // minimal level for quadrature

	minLevelLocal = minQuadratureLevel;          // minimal level for quadrature
	
  // 1. simplify the wavelets
  for (unsigned int i=0; i<waveletList.sizeWaveletList; ++i){
    // check if elements can be replaced by father and replace weight of
    // children by 0
		if (waveletList.W[i].level > minLevelLocal){
			noe = waveletList.W[i].noElements;      // number of entries in the element-list of the wavelet
			for (unsigned int s0=0; s0<noe; ++s0){
				j = waveletList.W[i].element[s0];   // element under consideration
        // check if it is a fine level element
				if (elementTree.element[j].level == waveletList.W[i].level) {
          // check if this is the 0. child in the father element
					j = elementTree.element[j].father;
					if (waveletList.W[i].element[s0] == elementTree.element[j].son[0]){
            // check 1st child
						for (s1=0; (s1<noe) && (elementTree.element[j].son[1]!=waveletList.W[i].element[s1]); ++s1);
						if ((s1 < noe) && (fabs(2*waveletList.W[i].weight[s0*noPhi+1]-waveletList.W[i].weight[s0*noPhi+0]-waveletList.W[i].weight[s1*noPhi+1]) < eps)){
              // check 2nd child
							for (s2=0; (s2<noe) && (elementTree.element[j].son[2] != waveletList.W[i].element[s2]); ++s2);
							if ((s2 < noe) && (fabs(2*waveletList.W[i].weight[s1*noPhi+2]-waveletList.W[i].weight[s1*noPhi+1]-waveletList.W[i].weight[s2*noPhi+2]) < eps)){
                // check 3rd child
								for (s3=0; (s3<noe) && (elementTree.element[j].son[3] != waveletList.W[i].element[s3]); ++s3); 
								{
									if ((s3 < noe) && (fabs(2*waveletList.W[i].weight[s2*noPhi+3]-waveletList.W[i].weight[s2*noPhi+2]-waveletList.W[i].weight[s3*noPhi+3]) < eps) \
												 && (fabs(2*waveletList.W[i].weight[s3*noPhi+0]-waveletList.W[i].weight[s3*noPhi+3]-waveletList.W[i].weight[s0*noPhi+0]) < eps) \
												 && (fabs(4*waveletList.W[i].weight[s0*noPhi+2]-waveletList.W[i].weight[s0*noPhi+0]-waveletList.W[i].weight[s1*noPhi+1]-waveletList.W[i].weight[s2*noPhi+2]-waveletList.W[i].weight[s3*noPhi+3]) < eps)) {
                    // all children have same weight, can be replaced
										waveletList.W[i].element[s0] = j;
										waveletList.W[i].weight[s0*noPhi+0] = 2*waveletList.W[i].weight[s0*noPhi+0];
										waveletList.W[i].weight[s0*noPhi+1] = 2*waveletList.W[i].weight[s1*noPhi+1];
										waveletList.W[i].weight[s0*noPhi+2] = 2*waveletList.W[i].weight[s2*noPhi+2];
										waveletList.W[i].weight[s0*noPhi+3] = 2*waveletList.W[i].weight[s3*noPhi+3];
										memset(&waveletList.W[i].weight[s1*noPhi],0,4*sizeof(double));
										memset(&waveletList.W[i].weight[s2*noPhi],0,4*sizeof(double));
										memset(&waveletList.W[i].weight[s3*noPhi],0,4*sizeof(double));
									}
								}
							}
						}
					}
				}
			}
   
      // adjust number of elements in wavelet
			k = 0;
			while (k < waveletList.W[i].noElements) {
				if ((waveletList.W[i].weight[k*noPhi+0] == 0) && (waveletList.W[i].weight[k*noPhi+1] == 0) && (waveletList.W[i].weight[k*noPhi+2] == 0) && (waveletList.W[i].weight[k*noPhi+3] == 0)){
					waveletList.W[i].noElements--;   // remove element
					for (unsigned int l=k; l<waveletList.W[i].noElements; ++l) {
						waveletList.W[i].element[l] = waveletList.W[i].element[l+1];
						memcpy(&waveletList.W[i].weight[l*noPhi],&waveletList.W[i].weight[(l+1)*noPhi],4*sizeof(double));
					}
				}
				else ++k;
			}

			waveletList.W[i].element = (unsigned int*) realloc(waveletList.W[i].element,waveletList.W[i].noElements*sizeof(unsigned int));
			waveletList.W[i].weight  = (double*) realloc(waveletList.W[i].weight, waveletList.W[i].noElements*4*sizeof(double));
		}
	}

  // 2. use prototypes for weights
	prototype = NULL;
	prototype_number = 0;

	for (unsigned int i=0; i<waveletList.sizeWaveletList; ++i){
		for (j=prototype_number-1; j>=0; j--){
			if (waveletList.W[i].noElements == waveletList.W[prototype[j]].noElements) {
        // same number of elements, check weights
				for (k=0; (k<waveletList.W[i].noElements) && (fabs(waveletList.W[i].weight[k*noPhi+0]-waveletList.W[prototype[j]].weight[k*noPhi+0]) < eps)
										 && (fabs(waveletList.W[i].weight[k*noPhi+1]-waveletList.W[prototype[j]].weight[k*noPhi+1]) < eps)
						 && (fabs(waveletList.W[i].weight[k*noPhi+2]-waveletList.W[prototype[j]].weight[k*noPhi+2]) < eps)
						 && (fabs(waveletList.W[i].weight[k*noPhi+3]-waveletList.W[prototype[j]].weight[k*noPhi+3]) < eps); ++k);

        // all weights are equal, replace weights with prototype
			  if (k == waveletList.W[i].noElements) {
				  free(waveletList.W[i].weight);
					waveletList.W[i].weight = waveletList.W[prototype[j]].weight;
					break;
				}
			}
		}

    // no prototype found, add to prototype list
		if (j == -1)	{
			if (prototype_number%delta == 0) prototype = (unsigned int*) realloc(prototype,(prototype_number+delta)*sizeof(unsigned int));
			prototype[prototype_number++] = i;
		}
	}

  // free memory for prototype indeces
	free(prototype);
	printf("%d prototypes\n",prototype_number);
#ifdef DEBUG
	FILE* debugFile = fopen("debug.out","a");
	fprintf(debugFile,">>> WAVELET_TREE_SIMPLIFY\n");
	for(unsigned int m = 0; m<waveletList.sizeWaveletList; ++m){
		fprintf(debugFile,"%d %d %d %d\n", m, waveletList.W[m].level, waveletList.W[m].noElements, waveletList.W[m].noSons);
		for(unsigned int i1 = 0 ; i1< waveletList.W[m].noElements; ++i1){
			fprintf(debugFile,"%d %lf %lf %lf %lf ", waveletList.W[m].element[i1], waveletList.W[m].weight[i1*4], waveletList.W[m].weight[i1*4+1], waveletList.W[m].weight[i1*4+2], waveletList.W[m].weight[i1*4+3]);
		}
		for(unsigned int i1 = 0 ; i1< waveletList.W[m].noSons; ++i1){
			fprintf(debugFile,"%d ", waveletList.W[m].son[i1]);
		}
		fprintf(debugFile,"\n");
	}
	fprintf(debugFile,"<<< WAVELET_TREE_SIMPLIFY\n");
	fclose(debugFile);
#endif

	return;
}

/// calculate if the interaction between wavelet ind1 and wavelet ind2 has to be
//computed
unsigned int LinAnsatzFunction::waveletWaveletCriterion(unsigned int ind1, unsigned int ind2, double c1, double c2){
  double dx, dy, dz;
  unsigned int i, j;
  unsigned int k,l;
  double h1, h2;
  double s2,t2, s1,t1;
  double dist;

#ifdef DEBUG
  FILE* debugFile = fopen("debug.out","a");
  fprintf(debugFile,"%d %d %lf %lf\n", ind1, ind2, c1, c2);
  fclose(debugFile);
#endif

  if(ind1 < ind2) return 0;

  // calculate distance between bounding boxes
  dx = fabs(B[ind1].mx - B[ind2].mx) - B[ind1].rx - B[ind2].rx;
  dy = fabs(B[ind1].my - B[ind2].my) - B[ind1].ry - B[ind2].ry;
  dz = fabs(B[ind1].mz - B[ind2].mz) - B[ind1].rz - B[ind2].rz;
  if (dx < 0) dx = 0;
  if (dy < 0) dy = 0;
  if (dz < 0) dz = 0;

  if( dx*dx+dy*dy+dz*dz >= c1*c1 ) return 0;

  // 2nd compression - test each element
  for(i = 0; i < waveletList.W[ind1].noElements;++i){
    k = waveletList.W[ind1].element[i];
  	h1 = 0.5/(1<<elementTree.element[k].level);
  	s1 = h1*(2*elementTree.element[k].index_s+1);
  	t1 = h1*(2*elementTree.element[k].index_t+1);
  	for(j = 0; j < waveletList.W[ind2].noElements; ++j){
  		l = waveletList.W[ind2].element[j];
  		if(elementTree.element[k].patch == elementTree.element[l].patch) {
  		  h2 = 0.5/(1 << elementTree.element[l].level);
  		  s2 = h2*(2*elementTree.element[l].index_s+1);
  		  t2 = h2*(2*elementTree.element[l].index_t+1);
  		  dist = fabs(s1-s2) < fabs(t1-t2) ? fabs(t1-t2) : fabs(s1-s2);
  		  if(fabs(fabs(dist-h2)-h1) < c2 ) return 1;
  		} 
      else {
  		  if ( distance(k, l) < c2) return 1;
  		}
  	}
  }
  return 0;
}

// memory release for bounding boxes
void LinAnsatzFunction::freeBoundingBoxes(){
  free(B);
  B = NULL;
}

// function that computes bounding boxes for wavelets
void LinAnsatzFunction::computeBoundingBoxes(){
  double minX, minY, minZ;
  double maxX, maxY, maxZ;

  unsigned int k;

  if(B) free(B);
  B = (BoundingBox*) malloc(waveletList.sizeWaveletList*sizeof(BoundingBox));
  for(unsigned int i = 0; i < waveletList.sizeWaveletList; ++i){
	  k = waveletList.W[i].element[0];
  	minX = elementTree.element[k].midpoint.x;
  	minY = elementTree.element[k].midpoint.y;
  	minZ = elementTree.element[k].midpoint.z;

  	maxX = elementTree.element[k].midpoint.x;
  	maxY = elementTree.element[k].midpoint.y;
  	maxZ = elementTree.element[k].midpoint.z;

  	for(unsigned int j = 0; j < waveletList.W[i].noElements; ++j){
	  	k = waveletList.W[i].element[j];
	  	if(elementTree.element[k].midpoint.x - elementTree.element[k].radius < minX) minX = elementTree.element[k].midpoint.x - elementTree.element[k].radius;
	  	if(elementTree.element[k].midpoint.y - elementTree.element[k].radius < minY) minY = elementTree.element[k].midpoint.y - elementTree.element[k].radius;
	  	if(elementTree.element[k].midpoint.z - elementTree.element[k].radius < minZ) minZ = elementTree.element[k].midpoint.z - elementTree.element[k].radius;

	  	if(elementTree.element[k].midpoint.x + elementTree.element[k].radius > maxX) maxX = elementTree.element[k].midpoint.x + elementTree.element[k].radius;
	  	if(elementTree.element[k].midpoint.y + elementTree.element[k].radius > maxY) maxY = elementTree.element[k].midpoint.y + elementTree.element[k].radius;
	  	if(elementTree.element[k].midpoint.z + elementTree.element[k].radius > maxZ) maxZ = elementTree.element[k].midpoint.z + elementTree.element[k].radius;
	  }

  	B[i].rx = 0.5*(maxX-minX);
  	B[i].ry = 0.5*(maxY-minY);
  	B[i].rz = 0.5*(maxZ-minZ);

  	B[i].mx = 0.5*(maxX+minX);
  	B[i].my = 0.5*(maxY+minY);
  	B[i].mz = 0.5*(maxZ+minZ);
  }
}

LinAnsatzFunction::~LinAnsatzFunction(){
  free(nodeList);
  for(unsigned int i = 0; i < elementTree.totalSizeElementList; ++i){
    free(elementTree.element[i].wavelet);
  }
  free(elementTree.element);
  for(unsigned int i = 0; i < nFunctions; ++i) free(pElementList[i]);
  free(pElementList);
  
  // delete prototypes
  unsigned int  *prototype;
  unsigned int  prototype_number;
  int j = 0;
  prototype = NULL;
  prototype_number = 0;

  for (unsigned int i=0; i<waveletList.sizeWaveletList; ++i) {
    for (j=prototype_number-1; j>=0; j--) {
      if (waveletList.W[i].weight == waveletList.W[prototype[j]].weight) {
        break;
      }
    }

    // prototype not found
    if (j == -1) {
      if (prototype_number%delta == 0) prototype = (unsigned int*) realloc(prototype,(prototype_number+delta)*sizeof(unsigned int));
      prototype[prototype_number++] = i;
    }
  }
  // free prototype memory
  for(unsigned int i = 0; i < prototype_number; ++i){
    free(waveletList.W[prototype[i]].weight);
  }
  free(prototype);

  for(unsigned int i = 0; i < waveletList.sizeWaveletList; ++i){
    free(waveletList.W[i].element);
    waveletList.W[i].weight = NULL;
    free(waveletList.W[i].son);
  }
  free(waveletList.W);
  delete(interCoeff);
  free(G);
}
///@todo implement tau, kappa, Phi times Phi, include memset

std::ostream & LinAnsatzFunction::printAnsatzFunction(std::ostream & os) 
{
  os << "A priori compression" << std::endl;      
  os << " a  = " << a << std::endl;
  os << " d' = " << dp << std::endl;
  os << "A posteriori compression" << std::endl;
  os << " b  = " << b;
  return os;
}
