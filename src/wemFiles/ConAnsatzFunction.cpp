#include "ConAnsatzFunction.hpp"
#include "Transformations.hpp"
#include "Mask.hpp"
#include "Vector3.hpp"
#include "Vector2.hpp"

#include <cstdio>
#include "string.h"

ConAnsatzFunction :: ConAnsatzFunction(){
  nLevels = 0;
  nFunctions = 0;
  nPatches = 0;
  minLevel = 1;
  noPhi = 1;

  B = NULL;
  B2 = NULL;

  td = 3;
  dp = 1.25;
  a = 1.25; ///< compression constant,  a > 1
  b = 0.001; ///< compression constant, 0 < b < 1

  quadratureLevel_ = 1;
  G = (SparseMatrix*) malloc(sizeof(SparseMatrix));
}

// NOTE the parameters must be initialized in the same order as declared in the header, otherwise a warning appears
ConAnsatzFunction :: ConAnsatzFunction(unsigned int _p, unsigned int _m, unsigned int _nf, Vector3*** pPointsIn){
  nLevels = _m;
  nFunctions = _nf;
  nPatches = _p;

  interCoeff  = new Interpolation(pPointsIn, 1, NEWTON, nLevels, nPatches);
  minLevel = 1;
  noPhi = 1;

  B = NULL;
  B2 = NULL;

  td = 3;
  dp = 1.25;
  a = 1.25; ///< compression constant,  a > 1
  b = 0.001; ///< compression constant, 0 < b < 1

  quadratureLevel_ = 1;
  G = (SparseMatrix*) malloc(sizeof(SparseMatrix));
}

ConAnsatzFunction :: ConAnsatzFunction(unsigned int _p, unsigned int _m, unsigned int _nf, double _a, double _b, double _dp, Vector3*** pPointsIn){
  nLevels = _m;
  nFunctions = _nf;
  nPatches = _p;

  interCoeff  = new Interpolation(pPointsIn, 1, NEWTON, nLevels, nPatches);
  minLevel = 1;
  noPhi = 1;

  B = NULL;
  B2 = NULL;

  td = 3;
  dp = _dp;
  a = _a; ///< compression constant,  a > 1
  b = _b; ///< compression constant, 0 < b < 1

  quadratureLevel_ = 1;
  G = (SparseMatrix*) malloc(sizeof(SparseMatrix));
}

void ConAnsatzFunction::setQuadratureLevel() {
  unsigned int  i, j, k, l;             // run index for wavelet/element list
  unsigned int  ind;                    // index of element under consideration
  unsigned int  minLevelLocal;          // miminal quadrature level
  unsigned int  noe;                    // number of elements of wavelet

  minLevelLocal = minQuadratureLevel;   // minimal quadrature level, global variable
  if (minLevelLocal > nLevels) minLevelLocal = nLevels;// check that the minimum level is less equal to the refinement level

  for (i=0; (i<waveletList.sizeWaveletList)&&(waveletList.W[i].level<=minLevelLocal); ++i) {
    for (j=waveletList.W[i].level; (j<=minLevelLocal); ++j) {
      noe = waveletList.W[i].noElements;
      for (k=0; k<noe; ++k) {
        // element under consideration
        ind = waveletList.W[i].element[k];
        // if the element is on a coarser level - refine 
        if (elementTree.element[ind].level < minLevelLocal) {
          waveletList.W[i].element[k] = elementTree.element[ind].son[0];
          waveletList.W[i].weight[k] = 0.5*waveletList.W[i].weight[k];
          waveletList.W[i].element = (unsigned int*) realloc(waveletList.W[i].element,(waveletList.W[i].noElements+3)*sizeof(unsigned int));
          waveletList.W[i].weight  = (double*) realloc(waveletList.W[i].weight ,(waveletList.W[i].noElements+3)*sizeof(double));

          // calculate weights
          for (l=1; l<4; ++l) {
            waveletList.W[i].element[waveletList.W[i].noElements] = elementTree.element[ind].son[l];
            waveletList.W[i].weight[waveletList.W[i].noElements] = waveletList.W[i].weight[k] ;
            waveletList.W[i].noElements++;
          }
        }
      }
    }
  }
#ifdef DEBUG
  FILE* debugFile = fopen("debug.out","a");
  fprintf(debugFile,">>> WAVELET_TREE_QUADRATURE\n");
  for(unsigned int m = 0; m<waveletList.sizeWaveletList; ++m){
    fprintf(debugFile,"%d %d %d\n", waveletList.W[m].level, waveletList.W[m].noElements, waveletList.W[m].noSons);
    for(unsigned int i1 = 0 ; i1< waveletList.W[m].noElements; ++i1){
      fprintf(debugFile,"%lf %lf ", waveletList.W[m].weight[i1], nodeList[elementTree.element[waveletList.W[m].element[i1]].vertex[0]].x);
    }
    fprintf(debugFile,"\n");
  }
  fprintf(debugFile,"<<< WAVELET_TREE_QUADRATURE\n");
  fflush(debugFile);
  fclose(debugFile);
#endif
  return;
}

void ConAnsatzFunction :: generateWaveletList(){
  unsigned int n = 1 << nLevels;  // nPatches*n*n elements on level m
  unsigned int  zw;                     // zero wavelet of i1 on level m
  unsigned int  ze;                     // zero element of i1 on level m
  SparseMatrix  Th, L;                  // mask matrices
  Wavelet     *w;                       // wavelet under consideration
  double newWeight;                     // new weights to be added
  const unsigned int sizeNewWeights = 1;// the size of the weights vector

  // length of the waveletList
  waveletList.sizeWaveletList = nPatches*n*n;
  //std::cout << "nw " << nw << std::endl;
  waveletList.W = (Wavelet*) calloc(waveletList.sizeWaveletList,sizeof(Wavelet));
  
  // loop over all levels till minimum
  for (unsigned int m=nLevels; m>=(signed int)minLevel; --m) {
    // compute wavelet masks
    dwtMask(&Th,&L,m,m,this);
    n = 1 << (m-1); // p*n*n elements on level m-1

    // 1. build scaling functions and ground wavelets?
    for (unsigned int i1=0; i1<nPatches; ++i1) {
      zw = i1*n*n;                    // zero wavelet of patch i1 on level m-1
      ze = nPatches*(4*n*n-1)/3+4*zw;  // zero element of patch i1 on level m

      for (unsigned int i2=0; i2<n; ++i2) {
        for (unsigned int i3=0; i3<n; ++i3) {  
          // 1. psi(s) * phi(t)
          w = &(waveletList.W[zw+nPatches*n*n+n*i2+i3]);
          w->level = m;               // level of wavelet
          w->noSons = 4;
          w->son = (unsigned int*) malloc(w->noSons*sizeof(unsigned int));
          w->son[0] = 4*nPatches*n*n+4*zw+4*i2*n+2*i3; // calculate children
          w->son[1] = w->son[0] + 1;
          w->son[2] = w->son[0] + 2*n;
          w->son[3] = w->son[0] + 2*n + 1;
          for (unsigned int s=0; s<Th.row_number[i3]; ++s) { // add scaling functions of level m
            newWeight = 0.5*Th.value[i3][s];
            addElement(w,elementTree,&newWeight,sizeNewWeights,ze+(4*n*i2  )+Th.index[i3][s]);
            addElement(w,elementTree,&newWeight, sizeNewWeights,ze+(4*n*i2+2*n)+Th.index[i3][s]);
          }

          // 2. phi_1(s) * psi(t) 
          w = &waveletList.W[zw+2*nPatches*n*n+n*i2+i3];
          w->level = m;               // level of wavelet
          w->noSons = 4;
          w->son = (unsigned int*) malloc(w->noSons*sizeof(unsigned int));
          w->son[0] = 8*nPatches*n*n+4*zw+4*i2*n+2*i3; // calculate children
          w->son[1] = w->son[0] + 2*n;
          w->son[2] = w->son[0] + 4*nPatches*n*n;
          w->son[3] = w->son[0] + 4*nPatches*n*n + 2*n;
          for (unsigned int t=0; t<Th.row_number[i2]; ++t){
            newWeight = 0.5*Th.value[i2][t];
            addElement(w,elementTree,&newWeight,sizeNewWeights,ze+(2*n*Th.index[i2][t])+(2*i3));
          }

          // 3. phi_2(s) * psi(t)
          w = &waveletList.W[zw+3*nPatches*n*n+n*i2+i3];
          w->level = m;               // level of wavelet
          w->noSons = 4;
          w->son = (unsigned int*) malloc(w->noSons*sizeof(unsigned int));
          w->son[0] = 8*nPatches*n*n+4*zw+4*i2*n+2*i3+1; // calculate children
          w->son[1] = w->son[0] + 2*n;
          w->son[2] = w->son[0] + 4*nPatches*n*n;
          w->son[3] = w->son[0] + 4*nPatches*n*n + 2*n;
          for (unsigned int t=0; t<Th.row_number[i2]; ++t){
            newWeight = 0.5*Th.value[i2][t];
            addElement(w, elementTree,&newWeight,sizeNewWeights,ze+(2*n*Th.index[i2][t])+(2*i3+1));
          }
        }
      }
    }
    freeSparse(&Th);
  }

  // add scaling functions of minLevel-1
  zw = 0;
  n = 1 << (minLevel-1);
  for (unsigned int i1=0; i1<nPatches; ++i1){
    for (unsigned int i2=0; i2<n; ++i2){
      for (unsigned int i3=0; i3<n; ++i3){
        waveletList.W[zw].level = minLevel-1;   // level for scaling function
        waveletList.W[zw].noSons = 4;
        waveletList.W[zw].son = (unsigned int*) malloc(waveletList.W[zw].noSons*sizeof(unsigned int));
        waveletList.W[zw].son[0] = zw;          // children of scaling function
        waveletList.W[zw].son[1] = zw+nPatches*n*n;
        waveletList.W[zw].son[2] = zw+2*nPatches*n*n;
        waveletList.W[zw].son[3] = zw+3*nPatches*n*n;
        newWeight = 1;
        addElement(&(waveletList.W)[zw],elementTree,&newWeight,sizeNewWeights,zw);   // element + enclosure
        ++zw;
      }
    }
  }
#ifdef DEBUG
  FILE* debugFile = fopen("debug.out","a");
  fprintf(debugFile,">>> WAVELET_TREE\n");
  for(unsigned int i2 = 0; i2<waveletList.sizeWaveletList; ++i2){
    fprintf(debugFile,"%d %d %d\n", waveletList.W[i2].level, waveletList.W[i2].noElements, waveletList.W[i2].noSons);
    for(unsigned int i1 = 0 ; i1< waveletList.W[i2].noElements; i1++){
      fprintf(debugFile,"%lf ", waveletList.W[i2].weight[i1]);
    }
    fprintf(debugFile,"\n");
  }
  fprintf(debugFile,"<<< WAVELET_TREE\n");
  fflush(debugFile);
  fclose(debugFile);
#endif
}

void ConAnsatzFunction::simplifyWaveletList(){
  unsigned int  k;               // run index
  signed int    j;               // run index
  unsigned int  ind;             // index wavelet under consideration
  unsigned int  s1, s2, s3;      // indeces of children
  unsigned int  noe;             // number of elements of wavelet
  unsigned int  *prototype;      // list for wavelet weight-prototypes
  unsigned int  prototype_number;// number of prototypes found
  unsigned int  minLevelLocal;   // minimal quadrature level, here copied from Constants

  minLevelLocal = minQuadratureLevel;

  // 1. simplify the wavelets
  for (unsigned int i=0; i<waveletList.sizeWaveletList; ++i){
    // check if elements can be replaced by father and replace weight of
    // children by 0
    if (waveletList.W[i].level > minLevelLocal) { 
      noe = waveletList.W[i].noElements;      // number of entries in the element-list of the wavelet
      for (unsigned int s0=0; s0<noe; ++s0) {
        ind = waveletList.W[i].element[s0];   // element under consideration
        // check if it is a fine level element
        if (elementTree.element[ind].level == waveletList.W[i].level) {
          // check if this is the 0. child in the father element
          ind = elementTree.element[ind].father;
          if ((int)waveletList.W[i].element[s0] == elementTree.element[ind].son[0]) {
            // check 1st child
            for (s1=0; (s1<noe) && (elementTree.element[ind].son[1]!=(int)waveletList.W[i].element[s1]); ++s1);
            if ((s1 < noe) && (waveletList.W[i].weight[s1] == waveletList.W[i].weight[s0])) {
              // check 2nd child
              for (s2=0; (s2<noe) && (elementTree.element[ind].son[2]!=(int)waveletList.W[i].element[s2]); ++s2);
              if ((s2 < noe) && (waveletList.W[i].weight[s2] == waveletList.W[i].weight[s0])) {
                // check 3rd child
                for (s3=0; (s3<noe) && (elementTree.element[ind].son[3]!=(int)waveletList.W[i].element[s3]); ++s3);
                if ((s3 < noe) && (waveletList.W[i].weight[s3] == waveletList.W[i].weight[s0])) {
                  // all children have same weight, can be replaced
                  waveletList.W[i].element[s0] = ind;
                  waveletList.W[i].weight[s0] *= 2;
                  waveletList.W[i].weight[s1] = waveletList.W[i].weight[s2] = waveletList.W[i].weight[s3] = 0;
                }
              }
            }
          }
        }
      }

      // adjust number of elements in wavelet
      k = 0;
      while (k < waveletList.W[i].noElements) {
        if (waveletList.W[i].weight[k] == 0) {
          waveletList.W[i].noElements--;   // remove element
          for (unsigned int l=k; l<waveletList.W[i].noElements; ++l) {
            waveletList.W[i].element[l] = waveletList.W[i].element[l+1];
            waveletList.W[i].weight[l] = waveletList.W[i].weight[l+1];
          }
        }
        else ++k;
      }
      waveletList.W[i].element = (unsigned int*) realloc(waveletList.W[i].element,waveletList.W[i].noElements*sizeof(unsigned int));
      waveletList.W[i].weight  = (double*)     realloc(waveletList.W[i].weight, waveletList.W[i].noElements*sizeof(double));
    }
  }

  // 2. use prototypes for weights
  prototype = NULL;
  prototype_number = 0;

  for (unsigned int i=0; i<waveletList.sizeWaveletList; ++i) {
    for (j=prototype_number-1; j>=0; j--) {
      if (waveletList.W[i].noElements == waveletList.W[prototype[j]].noElements) {
        // same number of elements, check weights
        for (k=0; (k<waveletList.W[i].noElements) && (waveletList.W[i].weight[k] == waveletList.W[prototype[j]].weight[k]); ++k);

        // all weights are equal, replace weights with prototype
        if (k == waveletList.W[i].noElements) {
          free(waveletList.W[i].weight);
          waveletList.W[i].weight = waveletList.W[prototype[j]].weight;
          break;
        }
      }
    }

    // no prototype found, add to prototype list
    if (j == -1) {
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
    fprintf(debugFile,"%d %d %d\n", waveletList.W[m].level, waveletList.W[m].noElements, waveletList.W[m].noSons);
    for(unsigned int i1 = 0 ; i1< waveletList.W[m].noElements; ++i1){
      fprintf(debugFile,"%lf %lf ", waveletList.W[m].weight[i1], nodeList[elementTree.element[waveletList.W[m].element[i1]].vertex[0]].x);
    }
    fprintf(debugFile,"\n");
  }
  fprintf(debugFile,"<<< WAVELET_TREE_SIMPLIFY\n");
  fflush(debugFile);
  fclose(debugFile);
#endif

  return;
}

/// calculate if the interaction between wavelet ind1 and wavelet ind2 has to be
//computed
unsigned int ConAnsatzFunction::waveletWaveletCriterion(unsigned int ind1, unsigned int ind2, double c1, double c2){
  double dx, dy, dz;
  unsigned int i, j;
  unsigned int k,l;
  double h1, h2;
  double s2, t2, s1, t1;
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

  // same patch, use in ConAF case the circles in the s-t plane
  if( B2[ind1].patch == B2[ind2].patch ){

    dx = fabs(B2[ind1].mx-B2[ind2].mx)-B2[ind1].rx-B2[ind2].rx;
    dy = fabs(B2[ind1].my-B2[ind2].my)-B2[ind1].ry-B2[ind2].ry;

    if( dx < 0 ) dx = 0;
    if( dy < 0 ) dy = 0;

    if( dx+dy > 0 ){
      if (dx*dx + dy*dy >= c2*c2 ) return 0;
      else return 1;
      //here in old code some for loop... why?
    } else {
      if((B2[ind1].mx - B2[ind1].rx - c2 < B2[ind2].mx - B2[ind2].rx) || \
        (B2[ind1].mx + B2[ind1].rx + c2 > B2[ind2].mx + B2[ind2].rx) || \
        (B2[ind1].my - B2[ind1].ry - c2 < B2[ind2].my - B2[ind2].ry) || \
        (B2[ind1].my + B2[ind1].ry + c2 > B2[ind2].my + B2[ind2].ry) )
        return 1;

      for( i = 0; i < waveletList.W[ind2].noElements; ++i){
        l = waveletList.W[ind2].element[i];
        h2 = 0.5/(1<<elementTree.element[l].level);
        s2 = h2*(2*elementTree.element[l].index_s+1);
        t2 = h2*(2*elementTree.element[l].index_t+1);
        dx = fabs(B2[ind1].mx - s2) - B2[ind1].rx - h2;
        dy = fabs(B2[ind1].my - t2) - B2[ind1].ry - h2;
        if (dx < 0) dx = 0;
        if (dy < 0) dy = 0;
        if( dx + dy > 0) {
          if (dx*dx+dy*dy < c2*c2) return 1;
        } else {
          if((B2[ind1].mx - B2[ind1].rx - c2 < s2-h2) || \
            (B2[ind1].mx + B2[ind1].rx + c2 > s2+h2) || \
            (B2[ind1].my - B2[ind1].ry - c2 < t2-h2) || \
            (B2[ind1].my + B2[ind1].ry + c2 > t2+h2)) {
            return 1;
          }
        }
      }
    }
    return 0;
  }

  // wavelets not on the same patch
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
      } else {
        if ( distance(k, l) < c2) return 1;
      }
    }
  }

  return 0;
}

// free method for the bounding boxes
void ConAnsatzFunction::freeBoundingBoxes(){
  free(B);
  free(B2);
  B = NULL;
  B2 = NULL;
}

// function that computes bounding boxes for wavelets
void ConAnsatzFunction::computeBoundingBoxes(){

  double minX, minY, minZ;
  double maxX, maxY, maxZ;

  unsigned int k;
  unsigned int l = nPatches;
  if(B) free(B);
  if(B2) free(B2);
  B  = (BoundingBox*) malloc(waveletList.sizeWaveletList*sizeof(BoundingBox));
  B2 = (BoundingBoxSquare*) malloc(waveletList.sizeWaveletList*sizeof(BoundingBoxSquare));
  for(unsigned int i = 0; i < waveletList.sizeWaveletList; ++i){
    k = waveletList.W[i].element[0];

    minX = elementTree.element[k].midpoint.x - elementTree.element[k].radius;
    minY = elementTree.element[k].midpoint.y - elementTree.element[k].radius;
    minZ = elementTree.element[k].midpoint.z - elementTree.element[k].radius;

    maxX = elementTree.element[k].midpoint.x + elementTree.element[k].radius;
    maxY = elementTree.element[k].midpoint.y + elementTree.element[k].radius;
    maxZ = elementTree.element[k].midpoint.z + elementTree.element[k].radius;

    // from all the elements in the support of the wavelet, choose min as
    // minimum and max as maximum
    for(unsigned int j = 0; j < waveletList.W[i].noElements; ++j){
      k = waveletList.W[i].element[j];
      if(elementTree.element[k].midpoint.x - elementTree.element[k].radius <minX) minX = elementTree.element[k].midpoint.x - elementTree.element[k].radius;
      if(elementTree.element[k].midpoint.y - elementTree.element[k].radius <minY) minY = elementTree.element[k].midpoint.y - elementTree.element[k].radius;
      if(elementTree.element[k].midpoint.z - elementTree.element[k].radius <minZ) minZ = elementTree.element[k].midpoint.z - elementTree.element[k].radius;

      if(elementTree.element[k].midpoint.x + elementTree.element[k].radius > maxX) maxX = elementTree.element[k].midpoint.x + elementTree.element[k].radius;
      if(elementTree.element[k].midpoint.y + elementTree.element[k].radius > maxY) maxY = elementTree.element[k].midpoint.y + elementTree.element[k].radius;
      if(elementTree.element[k].midpoint.z + elementTree.element[k].radius > maxZ) maxZ = elementTree.element[k].midpoint.z + elementTree.element[k].radius;
    }

    // construct bounding box radius in x, y, z direction
    B[i].rx = 0.5*(maxX-minX);
    B[i].ry = 0.5*(maxY-minY);
    B[i].rz = 0.5*(maxZ-minZ);

    // construct bounding box midpoint
    B[i].mx = 0.5*(maxX+minX);
    B[i].my = 0.5*(maxY+minY);
    B[i].mz = 0.5*(maxZ+minZ);
  }
#ifdef DEBUG
  FILE* debugFile = fopen("debug.out","a");
  fprintf(debugFile,">>> BBSQUARE\n");
  fclose(debugFile);
#endif

  // construct bounding boxes in s,t - space
  for(unsigned int i = 0; i < waveletList.sizeWaveletList; ++i){
    k = waveletList.W[i].element[0];
    B2[i].patch = elementTree.element[k].patch;

    minX = elementTree.element[k].index_s*1./(1<<elementTree.element[k].level);
    minY = elementTree.element[k].index_t*1./(1<<elementTree.element[k].level);

    maxX = (elementTree.element[k].index_s + 1)*1./(1<<elementTree.element[k].level);
    maxY = (elementTree.element[k].index_t + 1)*1./(1<<elementTree.element[k].level);

    for(unsigned int j = 1; j < waveletList.W[i].noElements; ++j){
      k = waveletList.W[i].element[j];

      if(elementTree.element[k].patch != B2[i].patch){
        B2[i].patch = l++;
        break;
      }

      if(elementTree.element[k].index_s < minX*(1<<elementTree.element[k].level)) minX = elementTree.element[k].index_s*1./(1<<elementTree.element[k].level);
      if(elementTree.element[k].index_t < minY*(1<<elementTree.element[k].level)) minY = elementTree.element[k].index_t*1./(1<<elementTree.element[k].level);

      if(elementTree.element[k].index_s + 1 > maxX*(1<<elementTree.element[k].level)) maxX = (elementTree.element[k].index_s+1)*1./(1<<elementTree.element[k].level);
      if(elementTree.element[k].index_t + 1 > maxY*(1<<elementTree.element[k].level)) maxY = (elementTree.element[k].index_t+1)*1./(1<<elementTree.element[k].level);
    }

    B2[i].rx = 0.5*(maxX-minX);
    B2[i].ry = 0.5*(maxY-minY);

    B2[i].mx = 0.5*(maxX+minX);
    B2[i].my = 0.5*(maxY+minY);
#ifdef DEBUG
  debugFile = fopen("debug.out","a");
  fprintf(debugFile,"%d %lf %lf %lf %lf\n", i, B2[i].rx, B2[i].ry, B2[i].mx, B2[i].my);
  fclose(debugFile);
#endif

  }
#ifdef DEBUG
  debugFile = fopen("debug.out","a");
  fprintf(debugFile,"<<< BBSQUARE\n");
  fclose(debugFile);
#endif

}

/**
 * helper function that stores into the specific vector of the AnsatzFunction
 * the values of the integral
 **/
void ConAnsatzFunction::calculateIntegral(double* a, double *c){
  c[0] = 0.5*(a[0] + a[3] + a[6] + a[ 9]);
  c[1] = 0.5*(a[1] + a[4] + a[7] + a[10]);
  c[2] = 0.5*(a[2] + a[5] + a[8] + a[11]);
}

/**
 * helper function that adds to the vector c the weight times the basis
 * functions evaluated at point xi
 */
void ConAnsatzFunction::calculateCRHS(double* c, double w, Vector2 xi){
  c[0] += w;
  return;
}

/**
 * helper function that stores into the specific vector of the AnsatzFunction
 * the contributions of the children
 */
void ConAnsatzFunction::calculateYRHS(double** y, int i){
  y[i][0] = 0.5*(y[elementTree.element[i].son[0]][0] + y[elementTree.element[i].son[1]][0] + y[elementTree.element[i].son[2]][0] + y[elementTree.element[i].son[3]][0]);
  return;
}

/**
 * helper function that stores into the energy vector the value at point xi
 **/
double ConAnsatzFunction::calculateUEnergy(double *u, Vector2 xi, unsigned int zi){
  return u[zi];
}

/**
 * function that creates the gramian matrix <phi_i, phi_j> - in this case, the
 * identity matrix
 */
void ConAnsatzFunction::createGram(unsigned int size, unsigned int maxRowNum){
  initSparse(G,size,size,maxRowNum);
  
  // create identyty mass matrix
  for (unsigned int i=0; i<size; ++i){
    addSparse(G,i,i,1);	
  }
}

// release the memory of the gramian
void ConAnsatzFunction::freeGram(){
  freeSparse(G);
}


/// inverse discrete wavelet transform of a
void ConAnsatzFunction :: tdwt(double *a){
  unsigned int n, zg;
  double*b;
  SparseMatrix T,L;

  b = (double*) calloc(waveletList.sizeWaveletList, sizeof(double));
  for(unsigned int m = minLevel; m <= nLevels; ++m){
    dwtMask(&T, &L, m, m, this);
    n = 1 << (m-1);

    // scaling functions and wavelets
    for(unsigned int i1 = 0; i1 < nPatches; ++i1){
      zg = i1*n*n;// zero element of patch i1
      for(unsigned int i2 = 0; i2 <n; ++i2){
        for(unsigned int i3 = 0; i3 < n; ++i3){
          b[4*zg+(4*n*i2    )+2*i3  ] += 0.5 * a[zg+n*i2+i3];
          b[4*zg+(4*n*i2    )+2*i3+1] += 0.5 * a[zg+n*i2+i3];
          b[4*zg+(4*n*i2+2*n)+2*i3  ] += 0.5 * a[zg+n*i2+i3];
          b[4*zg+(4*n*i2+2*n)+2*i3+1] += 0.5 * a[zg+n*i2+i3];
	    
          for (unsigned int s=0; s<T.row_number[i3]; ++s){
            b[4*zg+(4*n*i2    )+T.index[i3][s]] += 0.5 * T.value[i3][s] * a[zg+nPatches*n*n+n*i2+i3];
            b[4*zg+(4*n*i2+2*n)+T.index[i3][s]] += 0.5 * T.value[i3][s] * a[zg+nPatches*n*n+n*i2+i3];
          }
	       
          for (unsigned int t=0; t<T.row_number[i2]; ++t){
            b[4*zg+(2*n*T.index[i2][t])+(2*i3  )] += 0.5 * T.value[i2][t] * a[zg+2*nPatches*n*n+n*i2+i3];
            b[4*zg+(2*n*T.index[i2][t])+(2*i3+1)] += 0.5 * T.value[i2][t] * a[zg+3*nPatches*n*n+n*i2+i3];
          }
        }
      }
    }

    // 2. copy scaling functions in a
    memcpy(a,b,4*nPatches*n*n*sizeof(double));
    memset(b,0,4*nPatches*n*n*sizeof(double));
    freeSparse(&T);
  }

  // release memory of b
  free(b);
  return;
}

/// discrete wavelet transform of a
void ConAnsatzFunction :: dwt(double *a){
  unsigned int	n;      // n*n elements per patch
  unsigned int	zg;     // zero element on patch i1
  SparseMatrix		T, L; // mask matrices
  double		*b;         // help vector

  n = 1 << nLevels;
  b = (double*) calloc(waveletList.sizeWaveletList,sizeof(double));
  
  // loop over the levels
  for (unsigned int m=nLevels; m>=minLevel; --m){  
    dwtMask(&T,&L,m,nLevels,this);	// compute mask T and L
    n = 1 << (m-1);	// nPatches*n*n elements on level m-1
  
    // 1. scaling functions and wavelets
    for (unsigned int i1=0; i1<nPatches; ++i1){
      zg = i1*n*n;	// zero element on level m-1
      for (unsigned int i2=0; i2<n; ++i2){
        for (unsigned int i3=0; i3<n; ++i3){  
          // scaling functions
          b[zg+n*i2+i3] = 0.5 * ( a[4*zg+(4*n*i2    )+(2*i3  )] \
              + a[4*zg+(4*n*i2    )+(2*i3+1)] \
              + a[4*zg+(4*n*i2+2*n)+(2*i3  )] \
              + a[4*zg+(4*n*i2+2*n)+(2*i3+1)]);

          // wavelet: psi(s)*phi(t) 
          for (unsigned int s=0; s<T.row_number[i3]; ++s){
            b[zg+nPatches*n*n+n*i2+i3] += 0.5 * T.value[i3][s] * ( a[4*zg+(4*n*i2    )+T.index[i3][s]] \
                + a[4*zg+(4*n*i2+2*n)+T.index[i3][s]] );
          }

          // wavelet: phi_1(s)*psi(t) and phi_2(s)*psi(t)
          for (unsigned int t=0; t<T.row_number[i2]; ++t){
            b[zg+2*nPatches*n*n+n*i2+i3] += 0.5 * T.value[i2][t] * a[4*zg+(2*n*T.index[i2][t])+(2*i3  )];
            b[zg+3*nPatches*n*n+n*i2+i3] += 0.5 * T.value[i2][t] * a[4*zg+(2*n*T.index[i2][t])+(2*i3+1)];
          }
        }
      }
    }

    // 2. copy scaling functions to a
    memcpy(a,b,nPatches*n*n*sizeof(double));
    memset(b,0,nPatches*n*n*sizeof(double));
    freeSparse(&T);
  }

  memcpy(&a[nPatches*n*n],&b[nPatches*n*n],(waveletList.sizeWaveletList-nPatches*n*n)*sizeof(double));
  // free auxiliary memory
  free(b);
  return;
}

/**
 * helper function permutating the integral according to the AnsatzFunction used
 **/
void ConAnsatzFunction::permutate(double *a, double *b){
  a[0] = b[0];
  a[1] = b[2];
  a[2] = b[1];
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
void ConAnsatzFunction::quadratureGrade(signed int *g1, signed int*g2, int level1, int level2, double dist, double alpha) {
  alpha += level1+level2;
  dist = (dist < 1) ? log(dist) / log(2) : 0;

  // quadrature grade g1
  *g1 = (signed int) (-0.5*(alpha +dist)/(dist+level1+2));
  if(*g1 < 0) *g1 = 0;

  // quadrature grade g2
  *g2 = (signed int) (-0.5*(alpha +dist)/(dist+level2+2));
  if(*g1 < 0) *g2 = 0;

  return;
}

/// integration function for - no problem quadrature - modified scalar product
void ConAnsatzFunction::integrateNoProblem(double *c, unsigned int i1, unsigned int i2, unsigned int g1, Cubature *Q1, Cubature *Q2, double SingleLayer(Vector3, Vector3), double DoubleLayer(Vector3, Vector3, Vector3)) {
 	unsigned int i, j;
	double hs;
	Vector2 xi, eta;
	Vector3 y, n_y;
  double ht = 1./(1<<elementTree.element[i2].level);

  //calculate border integrals
  if(pRandWerte[g1].noP == 0){
    pRandWerte[g1].noP = Q1->noP;
    hs = 1./(1<<elementTree.element[i1].level);
    for(i = 0; i < Q1->noP; ++i){
      xi.x = hs*(elementTree.element[i1].index_s+Q1->xi[i].x);
      xi.y = hs*(elementTree.element[i1].index_t+Q1->xi[i].y);

      pRandWerte[g1].Chi[i] = interCoeff->Chi(xi, elementTree.element[i1].patch);
      //printf("intkon1 RW %d %lf %lf %lf %lf %lf %d %d %lf %lf %lf\n",i,hs, xi.x, xi.y, Q1->xi[i].x, Q1->xi[i].y, elementTree.element[i1].index_s, elementTree.element[i1].index_t, pRandWerte[g1].Chi[i].x, pRandWerte[g1].Chi[i].y,pRandWerte[g1].Chi[i].z); 
      pRandWerte[g1].n_Chi[i] = interCoeff->n_Chi(xi, elementTree.element[i1].patch);
      pRandWerte[g1].det_dChi[i] = hs * Q1->weight[i];
    }
  }

  // quadrature
  c[0] = c[1] = c[2] = 0;
  for (i = 0; i < Q2->noP; ++i) {
    eta.x = ht*(elementTree.element[i2].index_s + Q2->xi[i].x);
    eta.y = ht*(elementTree.element[i2].index_t + Q2->xi[i].y);
    y = interCoeff->Chi(eta,elementTree.element[i2].patch);
    n_y = interCoeff->n_Chi(eta,elementTree.element[i2].patch);
    for (j = 0; j < pRandWerte[g1].noP; ++j) {
      c[0] += Q2->weight[i] * pRandWerte[g1].det_dChi[j]*SingleLayer(pRandWerte[g1].Chi[j], y);
      c[1] += Q2->weight[i] * pRandWerte[g1].det_dChi[j]*DoubleLayer(pRandWerte[g1].Chi[j], y, n_y);
      c[2] += Q2->weight[i] * pRandWerte[g1].det_dChi[j]*DoubleLayer(y, pRandWerte[g1].Chi[j], pRandWerte[g1].n_Chi[j]);
    }
  }
  // L^2-normalized
  c[0] *= ht;
  c[1] *= ht;
  c[2] *= ht;
  //printf("intkon1 %lf %lf %lf\n",c[0], c[1], c[2]);
  return;
}

/// integration function for - same patches
void ConAnsatzFunction::integratePatch(double *c, unsigned int i1, Cubature * Q, double SingleLayer(Vector3, Vector3), double DoubleLayer(Vector3, Vector3, Vector3), double Identity) {
  unsigned int i, j;
  double w, d1, d2;
  double t1, t2, t3, t4;
  Vector2 eta, xi, a, b, s;
  Vector3 x, y;
  double h = 1./(1<<elementTree.element[i1].level);

  c[0] = c[1] = 0;
  s = Vector2 (h*elementTree.element[i1].index_s, h*elementTree.element[i1].index_t);

  for (i = 0; i < Q->noP; ++i) {
    xi = Q->xi[i];
    w = h * h * Q->weight[i] * xi.x * (1 - xi.x) * (1 - xi.x * xi.y);
    for (j = 0; j < Q->noP; ++j) {
      eta = Q->xi[j];
      t1 = h * eta.x * (1 - xi.x);
      t2 = h * eta.y * (1 - xi.x * xi.y);
      t3 = t1 + h * xi.x;
      t4 = t2 + h * xi.x * xi.y;

      //printf("INT PP00: %d %d\n", i, j);
      //printf("INT PP0: %lf %lf %lf %lf %lf\n", eta.x, eta.y, xi.x, xi.y, h);
      //printf("INT PP1: %lf %lf %lf %lf\n", t1, t2, t3, t4);
      
      a.x = s.x + t1;
      a.y = s.y + t2;
      b.x = s.x + t3;
      b.y = s.y + t4;
      x = interCoeff->Chi(a, elementTree.element[i1].patch);
      //printf("INT PP2: %d %lf %lf %lf %lf %lf %d %d\n", i1, a.x, a.y, x.x, x.y, x.z, nLevels, 1<<(nLevels-1));
      y = interCoeff->Chi(b, elementTree.element[i1].patch);
      //printf("INT PP3: %d %lf %lf %lf %lf %lf\n", i1, b.x, b.y, y.x, y.y, y.z);
      d1 = SingleLayer(x, y);
      d2 = DoubleLayer(x, y, interCoeff->n_Chi(b, elementTree.element[i1].patch))
          + DoubleLayer(y, x, interCoeff->n_Chi(a, elementTree.element[i1].patch));
      //printf("INT PPLAST1: %d %lf %lf\n", i1, d1, d2);

      a.y = s.y + t4;
      b.y = s.y + t2;
      x = interCoeff->Chi(a, elementTree.element[i1].patch);
      //printf("INT PP4: %d %lf %lf %lf %lf %lf\n", i1, a.x, a.y, x.x, x.y, x.z);
      y = interCoeff->Chi(b, elementTree.element[i1].patch);
      //printf("INT PP5: %d %lf %lf %lf %lf %lf\n", i1, b.x, b.y, y.x, y.y, y.z);
      d1 += SingleLayer(x, y);
      d2 += DoubleLayer(x, y, interCoeff->n_Chi(b, elementTree.element[i1].patch))
          + DoubleLayer(y, x, interCoeff->n_Chi(a, elementTree.element[i1].patch));
      //printf("INT PPLAST2: %d %lf %lf\n", i1, d1, d2);

      a.x = s.x + t2;
      a.y = s.y + t1;
      b.x = s.x + t4;
      b.y = s.y + t3;
      x = interCoeff->Chi(a, elementTree.element[i1].patch);
      //printf("INT PP6: %d %lf %lf %lf %lf %lf\n", i1, a.x, a.y, x.x, x.y, x.z);
      y = interCoeff->Chi(b, elementTree.element[i1].patch);
      //printf("INT PP7: %d %lf %lf %lf %lf %lf\n", i1, b.x, b.y, y.x, y.y, y.z);
      d1 += SingleLayer(x, y);
      d2 += DoubleLayer(x, y, interCoeff->n_Chi(b, elementTree.element[i1].patch))
          + DoubleLayer(y, x, interCoeff->n_Chi(a, elementTree.element[i1].patch));
      //printf("INT PPLAST3: %d %lf %lf\n", i1, d1, d2);

      a.y = s.y + t3;
      b.y = s.y + t1;
      x = interCoeff->Chi(a, elementTree.element[i1].patch);
      //printf("INT PP8: %d %lf %lf %lf %lf %lf\n", i1, a.x, a.y, x.x, x.y, x.z);
      y = interCoeff->Chi(b, elementTree.element[i1].patch);
      //printf("INT PP9: %d %lf %lf %lf %lf %lf\n", i1, b.x, b.y, y.x, y.y, y.z);
      d1 += SingleLayer(x, y);
      d2 += DoubleLayer(x, y, interCoeff->n_Chi(b, elementTree.element[i1].patch))
          + DoubleLayer(y, x, interCoeff->n_Chi(a, elementTree.element[i1].patch));
      //printf("INT PPLAST4: %d %lf %lf\n", i1, d1, d2);

      c[0] += 2 * w * Q->weight[j] * d1;
      c[1] += w * Q->weight[j] * d2;
    }
  }
  //printf("INT PPLAST: %d %lf %lf %lf %lf %lf %lf\n", i1, c[0], c[1], Identity, w, d1, d2);
  
  c[1] += Identity;
  c[2] = c[1];
  return;
}

/// integration function for - common edge
void ConAnsatzFunction::integrateEdge(double *c, unsigned int i1, unsigned int i2, unsigned int ind_s, unsigned int ind_t, Cubature *Q, double SingleLayer(Vector3, Vector3), double DoubleLayer(Vector3, Vector3, Vector3), double Identity) {
  unsigned int i, j;
  double w, d11, d12, d21, d22, d31, d32, t1, t2, t3, t4;
  Vector2 xi, eta, a, b, s, t;
  Vector3 x, y;
  double h = 1. / (1 << elementTree.element[i1].level);

  c[0] = c[1] = c[2] = 0;
  s = Vector2(h*elementTree.element[i1].index_s, h*elementTree.element[i1].index_t);
  t = Vector2(h*elementTree.element[i2].index_s, h*elementTree.element[i2].index_t);

  for (i = 0; i < Q->noP; ++i) {
    xi = Q->xi[i];
    w = xi.y * xi.y * Q->weight[i];
    t1 = xi.x * (1 - xi.y);
    t2 = (1 - xi.x) * (1 - xi.y);

    for (j = 0; j < Q->noP; ++j) {
      eta = vector2SMul(xi.y, Q->xi[j]);
      t3 = xi.x * (1 - eta.x);
      t4 = (1 - xi.x) * (1 - eta.x);

      a = kappa(s, tau(t1, eta.x, ind_s), h);
      b = kappa(t, tau(t2, eta.y, ind_t), h);
      x = interCoeff->Chi(a, elementTree.element[i1].patch);
      y = interCoeff->Chi(b, elementTree.element[i2].patch);
      d11 = SingleLayer(x, y);
      d21 = DoubleLayer(x, y, interCoeff->n_Chi(b, elementTree.element[i2].patch));
      d31 = DoubleLayer(y, x, interCoeff->n_Chi(a, elementTree.element[i1].patch));

      a = kappa(s, tau(1 - t1, eta.x, ind_s), h);
      b = kappa(t, tau(1 - t2, eta.y, ind_t), h);
      x = interCoeff->Chi(a, elementTree.element[i1].patch);
      y = interCoeff->Chi(b, elementTree.element[i2].patch);
      d11 += SingleLayer(x, y);
      d21 += DoubleLayer(x, y, interCoeff->n_Chi(b, elementTree.element[i2].patch));
      d31 += DoubleLayer(y, x, interCoeff->n_Chi(a, elementTree.element[i1].patch));

      a = kappa(s, tau(t3, xi.y, ind_s), h);
      b = kappa(t, tau(t4, eta.y, ind_t), h);
      x = interCoeff->Chi(a, elementTree.element[i1].patch);
      y = interCoeff->Chi(b, elementTree.element[i2].patch);
      d12 = SingleLayer(x, y);
      d22 = DoubleLayer(x, y, interCoeff->n_Chi(b, elementTree.element[i2].patch));
      d32 = DoubleLayer(y, x, interCoeff->n_Chi(a, elementTree.element[i1].patch));

      a = kappa(s, tau(1 - t3, xi.y, ind_s), h);
      b = kappa(t, tau(1 - t4, eta.y, ind_t), h);
      x = interCoeff->Chi(a, elementTree.element[i1].patch);
      y = interCoeff->Chi(b, elementTree.element[i2].patch);
      d12 += SingleLayer(x, y);
      d22 += DoubleLayer(x, y, interCoeff->n_Chi(b, elementTree.element[i2].patch));
      d32 += DoubleLayer(y, x, interCoeff->n_Chi(a, elementTree.element[i1].patch));

      a = kappa(s, tau(t4, eta.y, ind_s), h);
      b = kappa(t, tau(t3, xi.y, ind_t), h);
      x = interCoeff->Chi(a, elementTree.element[i1].patch);
      y = interCoeff->Chi(b, elementTree.element[i2].patch);
      d12 += SingleLayer(x, y);
      d22 += DoubleLayer(x, y, interCoeff->n_Chi(b, elementTree.element[i2].patch));
      d32 += DoubleLayer(y, x, interCoeff->n_Chi(a, elementTree.element[i1].patch));

      a = kappa(s, tau(1 - t4, eta.y, ind_s), h);
      b = kappa(t, tau(1 - t3, xi.y, ind_t), h);
      x = interCoeff->Chi(a, elementTree.element[i1].patch);
      y = interCoeff->Chi(b, elementTree.element[i2].patch);
      d12 += SingleLayer(x, y);
      d22 += DoubleLayer(x, y, interCoeff->n_Chi(b, elementTree.element[i2].patch));
      d32 += DoubleLayer(y, x, interCoeff->n_Chi(a, elementTree.element[i1].patch));

      c[0] += w * Q->weight[j] * ((1 - xi.y) * d11 + (1 - eta.x) * d12);
      c[1] += w * Q->weight[j] * ((1 - xi.y) * d21 + (1 - eta.x) * d22);
      c[2] += w * Q->weight[j] * ((1 - xi.y) * d31 + (1 - eta.x) * d32);
    }
  }

  // L^2-normalized
  c[0] *= h * h;
  c[1] *= h * h;
  c[2] *= h * h;
  return;
}

/// integration function for - common node in origin
void ConAnsatzFunction::integratePoint(double *c, unsigned int i1, unsigned int i2, unsigned int ind_s,
		unsigned int ind_t, Cubature *Q, double SingleLayer(Vector3, Vector3), double DoubleLayer(Vector3, Vector3, Vector3), double Identity) {
  unsigned int i, j;
  double w, d1, d2, d3;
  Vector2 xi, eta, a, s, t;
  Vector3 x1, n_x1, x2, n_x2, y1, n_y1, y2, n_y2, z, n_z;
  double h = 1. / (1 << elementTree.element[i1].level);

  c[0] = c[1] = c[2] = 0;
  s = Vector2(h*elementTree.element[i1].index_s, h*elementTree.element[i1].index_t);
  t = Vector2(h*elementTree.element[i2].index_s, h*elementTree.element[i2].index_t);

  for (i = 0; i < Q->noP; ++i) {
    xi = Q->xi[i];
    w = pow(xi.x, 3) * Q->weight[i];
    xi.y *= xi.x;

    a = kappa(s, tau(xi.x, xi.y, ind_s), h);
    x1 = interCoeff->Chi(a, elementTree.element[i1].patch);
    n_x1 = interCoeff->n_Chi(a, elementTree.element[i1].patch);

    a = kappa(s, tau(xi.y, xi.x, ind_s), h);
    x2 = interCoeff->Chi(a, elementTree.element[i1].patch);
    n_x2 = interCoeff->n_Chi(a, elementTree.element[i1].patch);

    a = kappa(t, tau(xi.x, xi.y, ind_t), h);
    y1 = interCoeff->Chi(a, elementTree.element[i2].patch);
    n_y1 = interCoeff->n_Chi(a, elementTree.element[i2].patch);

    a = kappa(t, tau(xi.y, xi.x, ind_t), h);
    y2 = interCoeff->Chi(a, elementTree.element[i2].patch);
    n_y2 = interCoeff->n_Chi(a, elementTree.element[i2].patch);

    for (j = 0; j < Q->noP; ++j) {
      eta = vector2SMul(xi.x, Q->xi[j]);

      a = kappa(t, tau(eta.x, eta.y, ind_t), h);
      z = interCoeff->Chi(a, elementTree.element[i2].patch);
      n_z = interCoeff->n_Chi(a, elementTree.element[i2].patch);

      d1 = SingleLayer(x1, z) + SingleLayer(x2, z);
      d2 = DoubleLayer(x1, z, n_z) + DoubleLayer(x2, z, n_z);
      d3 = DoubleLayer(z, x1, n_x1) + DoubleLayer(z, x2, n_x2);

      a = kappa(s, tau(eta.x, eta.y, ind_s), h);
      z = interCoeff->Chi(a, elementTree.element[i1].patch);
      n_z = interCoeff->n_Chi(a, elementTree.element[i1].patch);

      d1 += SingleLayer(z, y1) + SingleLayer(z, y2);
      d2 += DoubleLayer(z, y1, n_y1) + DoubleLayer(z, y2, n_y2);
      d3 += DoubleLayer(y1, z, n_z) + DoubleLayer(y2, z, n_z);

      c[0] += w * Q->weight[j] * d1;
      c[1] += w * Q->weight[j] * d2;
      c[2] += w * Q->weight[j] * d3;
    }
  }

  // L^2 - normalized
  c[0] *= h * h;
  c[1] *= h * h;
  c[2] *= h * h;
  return;
}

/// destructor - releases the memory assigned to this class
ConAnsatzFunction::~ConAnsatzFunction(){
  free(nodeList);
  for(unsigned int i = 0; i < elementTree.totalSizeElementList; ++i){
    free(elementTree.element[i].wavelet);
  }
  free(elementTree.element);
  for(unsigned int i = 0; i < nFunctions; ++i) {
    free(pElementList[i]);
  }
  free(pElementList);

  // delete prototypes
  unsigned int  *prototype;
  unsigned int  prototype_number;
  int j = 0;
  prototype = NULL;
  prototype_number = 0;

  for (unsigned int i=0; i<waveletList.sizeWaveletList; i++) {
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
