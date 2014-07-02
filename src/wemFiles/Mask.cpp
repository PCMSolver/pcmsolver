/**
 * @file Mask.cpp
 *
 * @brief defines the fast wavelet transform masks for constant with 3
 * vanishin moments and linear wavelets with 4 vanishing moments and
 * double points at patch edges
 */

#include "Mask.hpp"
//#include "LinAnsatzFunction.hpp"
//#include "ConAnsatzFunction.hpp"
class ConAnsatzFunction;
class LinAnsatzFunction;
 
// constant - defines the mask T on level m
void mask_T11(SparseMatrix * T, unsigned int m) {
  unsigned int    n = 1 << (m-1);
  unsigned int    i;

  initSparse(T,n,2*n,2);

  for (i=0; i<n; i++) {
    setSparse(T,i,2*i  ,+1);
    setSparse(T,i,2*i+1,-1);
  }
  return;
}

// constant - defines the mask T on level m
void mask_T13(SparseMatrix* T, unsigned int m) {
  unsigned int    n = 1 << (m-1);
  unsigned int    i;

  initSparse(T,n,2*n,6);

  for (i=0; i<n; i++) {
    if (i == 0) {
      // left boundary wavelet
      setSparse(T,0,0,+0.625);
      setSparse(T,0,1,-1.375);
      setSparse(T,0,2,+0.500);
      setSparse(T,0,3,+0.500);
      setSparse(T,0,4,-0.125);
      setSparse(T,0,5,-0.125);
    } else if (i < n-1) {
      // strationary wavelet
      setSparse(T,i,2*i-1,-0.25);
      setSparse(T,i,2*i  ,+0.75);
      setSparse(T,i,2*i+1,-0.75);
      setSparse(T,i,2*i+2,+0.25);
    } else if (i == n-1) {
      // right boundary wavelet
      setSparse(T,n-1,2*n-6,+0.125);
      setSparse(T,n-1,2*n-5,+0.125);
      setSparse(T,n-1,2*n-4,-0.500);
      setSparse(T,n-1,2*n-3,-0.500);
      setSparse(T,n-1,2*n-2,+1.375);
      setSparse(T,n-1,2*n-1,-0.625);
    }
  }
  return;
}


// linear - defines the mask T on level m, 2 vanishing moments
void mask_T22(SparseMatrix *T, unsigned int m, unsigned int M) {
  unsigned int    n = 1 << m;
  unsigned int    i;

  initSparse(T,n+1,n+1,4);

  for (i=0; i<=n; i++) {
    if (i%2 == 0) { 
      // scaling functions
      if (i == 0) {
        // left scaling function
        setSparse(T,0,0,+1.0);
        setSparse(T,0,1,+0.5);
      } else if (i < n) {
        // stationary scaling function
        setSparse(T,i,i-1,0.5);
        setSparse(T,i,i  ,1.0);
        setSparse(T,i,i+1,0.5);
      } else if (i == n) {
        // right scaling function
        setSparse(T,n,n-1,+0.5);
        setSparse(T,n,n  ,+1.0);
      }
    } else {
      // wavelet
      if (i == 1) {
        // left boundary wavelet
        setSparse(T,1,0,-0.7500);
        setSparse(T,1,1,+0.5625);
        setSparse(T,1,2,-0.1250);
        setSparse(T,1,3,-0.0625);
      } else if (i < n-1) {
        // stationary wavelet
        setSparse(T,i,i-1,-0.5);
        setSparse(T,i,i  ,+1.0);
        setSparse(T,i,i+1,-0.5);
      } else if (i == n-1) {
        // right boundary wavelet
        setSparse(T,n-1,n-3,-0.0625);
        setSparse(T,n-1,n-2,-0.1250);
        setSparse(T,n-1,n-1,+0.5625);
        setSparse(T,n-1,n  ,-0.7500);
      }
    }
  }
  return;
}

// linear - defines the mask T on level m, 4 vanishing moments
void mask_T24(SparseMatrix *T, unsigned int m, unsigned int M) {
  unsigned int    n = 1 << m;
  unsigned int    i;

  initSparse(T,n+1,n+1,8);

  for (i=0; i<=n; i++) {
    if (i%2 == 0) {
      // scaling functions
      if (i == 0) {
        // left scaling function
        setSparse(T,0,0,+1.0);
        setSparse(T,0,1,+0.5);
      } else if (i < n) {
        // stationary scaling function
        setSparse(T,i,i-1,0.5);
        setSparse(T,i,i  ,1.0);
        setSparse(T,i,i+1,0.5);
      } else if (i == n) {
        // right scaling function
        setSparse(T,n,n-1,+0.5);
        setSparse(T,n,n  ,+1.0);
      }
    } else {
      // wavelet
      if (i == 1) {
        // left boundary wavelet
        setSparse(T,1,0,-35./64);
        setSparse(T,1,1,875./1536);
        setSparse(T,1,2,-241./768);
        setSparse(T,1,3,-53./512);
        setSparse(T,1,4,41./384);
        setSparse(T,1,5,67./1536);
        setSparse(T,1,6,-5./256);
        setSparse(T,1,7,-5./512);
      } else if (i < n-1) {
        // stationary wavelet
        setSparse(T,i,i-2,+0.125);
        setSparse(T,i,i-1,-0.500);
        setSparse(T,i,i  ,+0.750);
        setSparse(T,i,i+1,-0.500);
        setSparse(T,i,i+2,+0.125);
      } else if (i == n-1) {
        // right boundary wavelet
        setSparse(T,n-1,n-7,-5./512);
        setSparse(T,n-1,n-6,-5./256);
        setSparse(T,n-1,n-5,67./1536);
        setSparse(T,n-1,n-4,41./384);
        setSparse(T,n-1,n-3,-53./512);
        setSparse(T,n-1,n-2,-241./768);
        setSparse(T,n-1,n-1,875./1536);
        setSparse(T,n-1,n  ,-35./64);
      }
    }
  }
  return;
}


// choose the correct mask according to level
template <> void dwtMask(SparseMatrix *T, SparseMatrix *L, unsigned int m, unsigned int M, LinAnsatzFunction *linAF) {
  minLevel = 2;
  if (m <= 2) mask_T22(T,m,M);
  else        mask_T24(T,m,M);
  return;
}

// choose the correct mask according to level
template <> void dwtMask(SparseMatrix *T, SparseMatrix *L, unsigned int m, unsigned int M, ConAnsatzFunction *conAF) {
  minLevel = 1;
  if (m < 3) {  mask_T11(T,m);  }
  else       {  mask_T13(T,m);  }
  return;
}

