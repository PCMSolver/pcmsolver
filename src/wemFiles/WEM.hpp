/**
 * @file WEM.hpp
 * computes the stffnessmatrix for the Galerkin method using a template parameter for the ansatz function
 */
 
#ifndef WEM_HPP
#define WEM_HPP

#if defined DEBUG || defined DEBUG2
#include <stdio.h>
#endif

#include <string.h>

#include "Cubature.hpp"
#include "GaussSquare.hpp"
#include "GenericAnsatzFunction.hpp"

class Vector3;

/**
 * @param af the ansatz function that is used at the moment, can be ConAnsatzFunction or LinAnsatzFunction
 * @param sM the sparse matrix that is computed
 * @param SingleLayer the single layer function to be used in this call
 * @param DoubleLayer the double layer function that is called here
 * @param Identity the normalization factor?
*/
template <class T> void WEM(T* af, SparseMatrix *sM, double SingleLayer(Vector3, Vector3), double DoubleLayer(Vector3, Vector3, Vector3), double Identity) {
  unsigned int *ansatzWavelet;
  unsigned int *weightIndex; // double *weight; //in ConAF
  unsigned int testWavelet;
  Cubature *Q;
  double **prec;
  double *c;
  double *s;
  double *t;

  unsigned int sizeT = 3, sizeS = sizeT*af->noPhi, sizeC = sizeS*af->noPhi;

  // Initialization
  //printf("WEM IN %d %d %d\n",af->sizeC, af->sizeS, af->sizeT);
  c = (double*) calloc(sizeC,sizeof(double));
  s = (double*) calloc(sizeS,sizeof(double));
  t = (double*) calloc(sizeT,sizeof(double));
  weightIndex = NULL;
  ansatzWavelet = NULL;

  af->initRandwerte(g_max);
  initGaussSquare(&Q,g_max); // Kubatur-Formeln
#if defined DEBUG || defined DEBUG2
  FILE *debugFile;
#endif
#ifdef DEBUG
  debugFile = fopen("debug.out","a");
  fprintf(debugFile, ">>> CUBATURE\n");
  for(int i = 0; i< g_max; i++) {
    fprintf(debugFile,"%d\n",i);
      for(unsigned j = 0; j < Q[i].noP;++j)
        fprintf(debugFile, "%d %f %f %f\n", i, Q[i].xi[j].x, Q[i].xi[j].y, Q[i].weight[j]);
  }
  fprintf(debugFile, "<<< CUBATURE\n");
  fclose(debugFile);
#endif
  
  // calculate log2 of precision
  prec = (double **) malloc((af->nLevels+1)*sizeof(double*) + (af->nLevels+1)*(af->nLevels+1)*sizeof(double));
  for(int i = 0; i <= (int)af->nLevels; ++i) {
    prec[i] = (double*)(prec+(af->nLevels+1))+i*(af->nLevels+1);
    for(int  j = 0; j <= (int)af->nLevels; ++j){
      prec[i][j] = (2*af->nLevels-(i+j))*(op-2*af->dp)/(2*af->td+op);
      if( prec[i][j] > -fabs(i-j) ) prec[i][j] = -fabs(i-j);
      prec[i][j] += op*af->nLevels-af->dp*(2*af->nLevels-(i+j));
    }
  }
#ifdef DEBUG
  debugFile = fopen("debug.out","a");
  fprintf(debugFile, ">>> PREC\n");
  for(unsigned int i = 0; i< af->nLevels+1;i++) {
    fprintf(debugFile,"%d\n",i);
      for(unsigned j = 0; j < af->nLevels+1;++j)
        fprintf(debugFile, "%d %lf\n", j, prec[i][j]);
  }
  fprintf(debugFile, "<<< PREC\n");
  fclose(debugFile);
#endif
  // calculate system matrix linewise
  for(unsigned int i = 0; i < af->totalSizeElementList; ++i){
    af->resetRandwerte(g_max);
    weightIndex = (unsigned int*) realloc(weightIndex, af->elementTree.element[i].noWavelets*sizeof(unsigned int));
    ansatzWavelet = (unsigned int*) realloc(ansatzWavelet, af->elementTree.element[i].noWavelets*sizeof(unsigned int));
    memset(ansatzWavelet,0, af->elementTree.element[i].noWavelets*sizeof(unsigned int));

    // determine first testWavelet and the weights of the ansatzFunction
    // regarding element i
    testWavelet = sM->n;
    for(unsigned int j = 0; j < af->elementTree.element[i].noWavelets; ++j){
      if(sM->index[af->elementTree.element[i].wavelet[j]][0] < testWavelet) testWavelet = sM->index[af->elementTree.element[i].wavelet[j]][0];
      for(weightIndex[j] = 0; af->waveletList.W[af->elementTree.element[i].wavelet[j]].element[weightIndex[j]] != i; weightIndex[j]++);
    }

    // loop over the lines of the System Matrix sM
    while((testWavelet < sM->n) && (af->elementTree.element[i].level+2 > af->waveletList.W[testWavelet].level)){
      // test quadrature of wavelet i with testWavelet
      memset(s, 0, sizeS*sizeof(double)); // sizeS = 12 LinAF 3 ConAF
      for(unsigned int j = 0; j < af->waveletList.W[testWavelet].noElements; ++j){
        if( i > af->waveletList.W[testWavelet].element[j] ) {
          af->elementElementInteraction(c, i, af->waveletList.W[testWavelet].element[j], prec[af->elementTree.element[i].level][af->waveletList.W[testWavelet].level], Q, SingleLayer, DoubleLayer, Identity); 
          for(unsigned int k = 0; k < sizeS; ++k){
            for(unsigned int k2 = 0; k2 < af->noPhi;++k2){
              s[k] += af->waveletList.W[testWavelet].weight[j*af->noPhi+k2] * c[k*af->noPhi+k2];
              //printf("ELELINT1a %d %d %d %d %d %lf %lf %lf\n", i, k, k2, testWavelet, af->waveletList.W[testWavelet].element[j], af->waveletList.W[testWavelet].weight[j*af->noPhi+k2],c[k*af->noPhi+k2], s[k]);
              //printf("ELELINT1a %d %d %d %lf %lf %lf\n", k, k2, af->noPhi, 0.0, 0.0, 0.0);
            }
          }
        } else if (i == af->waveletList.W[testWavelet].element[j]){
          af->elementElementInteraction(c, i, af->waveletList.W[testWavelet].element[j], prec[af->elementTree.element[i].level][af->waveletList.W[testWavelet].level], Q, SingleLayer, DoubleLayer, Identity); 
          for( unsigned int k = 0; k < sizeS; ++k){
            for(unsigned int k2 = 0; k2 < af->noPhi;++k2){
              s[k] += 0.5 * af->waveletList.W[testWavelet].weight[j*af->noPhi+k2] * c[k*af->noPhi+k2];
              //printf("ELELINT1b %d %d %d %d %d %lf %lf %lf %lf\n", i, k, k2, testWavelet, af->waveletList.W[testWavelet].element[j], af->waveletList.W[testWavelet].weight[j*af->noPhi+k2],c[k*af->noPhi+k2], s[k], Identity);
            }
          }
        }
      }
      // save calculated integrals
      for(unsigned int j = 0; j < af->elementTree.element[i].noWavelets; ++j){
        if( sM->index[af->elementTree.element[i].wavelet[j]][ansatzWavelet[j]] == testWavelet){
          memset(t, 0, sizeT*sizeof(double));
          for (unsigned int k = 0; k < sizeT; ++k){
            for(unsigned int k2 = 0; k2 < af->noPhi; ++k2){
              t[k] += af->waveletList.W[af->elementTree.element[i].wavelet[j]].weight[weightIndex[j]*af->noPhi+k2] * s[k*af->noPhi+k2];
              //printf("ELELINT2 %d %lf %lf %lf\n", k, af->waveletList.W[af->elementTree.element[i].wavelet[j]].weight[weightIndex[j]*af->noPhi+k2], s[k*af->noPhi+k2], t[k]);
            }
          }
          //printf("ELELINT3 %d %d %d %d %lf\n", i, j, af->elementTree.element[i].wavelet[j], ansatzWavelet[j], t[0]);
          sM->value1[af->elementTree.element[i].wavelet[j]][ansatzWavelet[j]] +=t[0];
          sM->value2[af->elementTree.element[i].wavelet[j]][ansatzWavelet[j]] +=t[1];
          addSparse2(sM, testWavelet, af->elementTree.element[i].wavelet[j], t[0], t[2]);
          ++ansatzWavelet[j];
        }
      }
      // determine new testWavelet
      testWavelet = sM->n;
      for(unsigned int j = 0; j < af->elementTree.element[i].noWavelets; ++j){
        if (sM->index[af->elementTree.element[i].wavelet[j]][ansatzWavelet[j]] < testWavelet)
          testWavelet = sM->index[af->elementTree.element[i].wavelet[j]][ansatzWavelet[j]];
      }
    }
    // erase element element interactions
    for(unsigned int j = 0; j < af->elementTree.element[i].interaction.integralNumber; ++j) free(af->elementTree.element[i].interaction.value[j]);
    af->elementTree.element[i].interaction.integralNumber = 0;
    free(af->elementTree.element[i].interaction.value);
    free(af->elementTree.element[i].interaction.index);
#ifdef DEBUG2
  debugFile = fopen("debug.out","a");
  fprintf(debugFile,">>> RANDWERTE %d\n",i);
  for(int i1 = 0; i1< g_max; ++i1) {
    fprintf(debugFile,"%d %d\n",i1, af->pRandWerte[i1].noP);
      for(unsigned int i2 = 0; i2 < af->pRandWerte[i1].noP; ++i2)
        fprintf(debugFile, "%d %d %f %f %f %f %f %f %f\n", i1, i2, af->pRandWerte[i1].Chi[i2].x, af->pRandWerte[i1].Chi[i2].y, af->pRandWerte[i1].Chi[i2].z, af->pRandWerte[i1].n_Chi[i2].x, af->pRandWerte[i1].n_Chi[i2].y, af->pRandWerte[i1].n_Chi[i2].z, af->pRandWerte[i1].det_dChi[i2]);
  }
  fprintf(debugFile,"<<< RANDWERTE %d\n",i);
  fclose(debugFile);
#endif
  }
  // release memory
  free(prec);
  freeGaussSquare(&Q, g_max);
  af->freeRandwerte();
  free(ansatzWavelet);
  free(weightIndex);
  free(s);
  free(t);
  free(c);
}
#endif
