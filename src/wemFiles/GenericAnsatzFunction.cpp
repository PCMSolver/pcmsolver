#include "GenericAnsatzFunction.hpp"
#include "Topology.hpp"
#include "Vector3.hpp"
#include "Vector2.hpp"

#include "stdio.h"
#include <string.h>
#include <math.h>

/**
 * @todo check if one can allocate memory in one go
 * @todo this return will be absolete
 */
unsigned int GenericAnsatzFunction :: genNet(Vector3 ***U) {
  unsigned int i;
  unsigned int n = 1<<nLevels; // n
  et_node *pF;

  double r1,r2;        // help values for radius of circumscribed sphere
  Vector3 d1, d2;      // help values for midpoint of circumscribed sphere

  unsigned int  nLevel;// n = 2^(m-1)
  unsigned int  S;     // distance to next point on current level
  unsigned int  fz;    // element counter

  // allocate initial memory
  nodeList = (Vector3*) malloc(4*nPatches*sizeof(Vector3));
  pElementList = (unsigned int**) malloc(nPatches*sizeof(unsigned int*));

  elementTree.nop = nPatches;
  elementTree.totalSizeElementList = nPatches*(4*(1<<2*nLevels)-1)/3;
  //allocate all the nodes at once
  pF = elementTree.element = (et_node*) calloc(elementTree.totalSizeElementList, sizeof(et_node));

  // calculate edge point of parameter domain
  nNodes = 4*nPatches;
  nFunctions = nPatches;
  for (i=0; i<nPatches; ++i) {
    pElementList[i] = (unsigned int*) calloc(4,sizeof(unsigned int));

    pF[i].level = 0;
    pF[i].number = i;
    pF[i].patch = i;
    pF[i].index_s = 0;
    pF[i].index_t = 0;
    pF[i].father = -1;

    pF[i].vertex[0] = 4*i;
    pF[i].vertex[1] = 4*i+1;
    pF[i].vertex[2] = 4*i+3;
    pF[i].vertex[3] = 4*i+2;

    pF[i].son[0] = nPatches+4*i;
    pF[i].son[1] = nPatches+4*i+1;
    pF[i].son[2] = nPatches+4*i+3;
    pF[i].son[3] = nPatches+4*i+2;

    pF[nPatches+4*i  ].father = i;
    pF[nPatches+4*i+1].father = i;
    pF[nPatches+4*i+2].father = i;
    pF[nPatches+4*i+3].father = i;

    // calculate 4 corner points
    nodeList[4*i]   = U[i][0][0];
    nodeList[4*i+1] = U[i][0][n]; 
    nodeList[4*i+2] = U[i][n][0]; 
    nodeList[4*i+3] = U[i][n][n]; 

    pElementList[i][0] = 4*i;
    pElementList[i][1] = 4*i+1;
    pElementList[i][2] = 4*i+3;
    pElementList[i][3] = 4*i+2;
  }

  for (unsigned int auxLevel=1; auxLevel<=nLevels; ++auxLevel){
    nLevel = 1 << (auxLevel-1);   // n = 2^(m-1)
    S = 1 << (nLevels-auxLevel);  // distance to next point on current level

    // allocate memory
    nodeList = (Vector3*) realloc(nodeList,4*(nNodes)*sizeof(Vector3));
    pElementList = (unsigned int**) realloc(pElementList,4*(nFunctions)*sizeof(unsigned int*));
    for (unsigned int i1=nFunctions; i1<4*(nFunctions); ++i1) pElementList[i1] = (unsigned int*) calloc(4,sizeof(unsigned int));

    // copy old elements into new level
    fz = (nPatches-1)*nLevel*nLevel;
    i = nPatches*(4*nLevel*nLevel-1)/3;
    for (int i1=nPatches-1; i1>=0; i1--) {
      for (int i2=nLevel-1; i2>=0; i2--) {
        for (int i3=nLevel-1; i3>=0; i3--) {
          pF[i+4*fz+2*nLevel*(2*i2+1)+2*i3  ].vertex[3] = pElementList[fz+nLevel*i2+i3][3];
          pElementList[4*fz+2*nLevel*(2*i2+1)+2*i3  ][3] = pElementList[fz+nLevel*i2+i3][3];

          pF[i+4*fz+2*nLevel*(2*i2+1)+2*i3+1].vertex[2] = pElementList[fz+nLevel*i2+i3][2];
          pElementList[4*fz+2*nLevel*(2*i2+1)+2*i3+1][2] = pElementList[fz+nLevel*i2+i3][2];

          pF[i+4*fz+4*nLevel*i2    +2*i3+1].vertex[1] = pElementList[fz+nLevel*i2+i3][1];
          pElementList[4*fz+4*nLevel*i2    +2*i3+1][1] = pElementList[fz+nLevel*i2+i3][1];

          pF[i+4*fz+4*nLevel*i2    +2*i3  ].vertex[0] = pElementList[fz+nLevel*i2+i3][0];
          pElementList[4*fz+4*nLevel*i2    +2*i3  ][0] = pElementList[fz+nLevel*i2+i3][0];

          pF[i+4*fz+4*nLevel*i2    +2*i3  ].level = auxLevel;
          pF[i+4*fz+2*nLevel*2*i2    +2*i3  ].patch = i1;
          pF[i+4*fz+2*nLevel*2*i2    +2*i3  ].father =  nPatches*(nLevel*nLevel-1)/3 + i1*nLevel*nLevel+ i2*nLevel + i3;

          pF[i+4*fz+2*nLevel*2*i2    +2*i3  ].son[0] =  nPatches*(4*nLevel*4*nLevel-1)/3 + 16*i1*nLevel*nLevel + 4*i2*4*nLevel + 4*i3;
          pF[i+4*fz+2*nLevel*2*i2    +2*i3  ].son[1] =  nPatches*(4*nLevel*4*nLevel-1)/3 + 16*i1*nLevel*nLevel + 4*i2*4*nLevel + 4*i3+1;
          pF[i+4*fz+2*nLevel*2*i2    +2*i3  ].son[2] =  nPatches*(4*nLevel*4*nLevel-1)/3 + 16*i1*nLevel*nLevel + 4*(4*i2+1)*nLevel + 4*i3+1;
          pF[i+4*fz+2*nLevel*2*i2    +2*i3  ].son[3] =  nPatches*(4*nLevel*4*nLevel-1)/3 + 16*i1*nLevel*nLevel + 4*(4*i2+1)*nLevel + 4*i3;

          pF[i+4*fz+4*nLevel*i2    +2*i3  ].number = pF[pF[i+4*fz+4*nLevel*i2    +2*i3  ].father].number*4;

          pF[i+4*fz+4*nLevel*i2    +2*i3  ].index_s = 2*i3;
          pF[i+4*fz+4*nLevel*i2    +2*i3  ].index_t = 2*i2;
        }
      }
      fz -= nLevel*nLevel;
    }

    // add new points and elements for current level
    fz = 0;
    nFunctions *= 4;
    for (unsigned int i1=0; i1<nPatches; ++i1) {
      for (unsigned int i2=0; i2<=2*nLevel; ++i2) {
        for (unsigned int i3=0; i3<=2*nLevel; ++i3) {

          if ((i2%2 == 1) || (i3%2 == 1)) {
            nodeList[nNodes] = U[i1][i2*S][i3*S];
            if ((i2<2*nLevel) && (i3<2*nLevel)){
              pElementList[fz+2*nLevel* i2   +i3  ][0] = nNodes;
              pF[i+ fz+2*nLevel* i2 +i3  ].vertex[0] = nNodes;

              pF[i+ fz+2*nLevel* i2 +i3  ].level = auxLevel;
              pF[i+ fz+2*nLevel* i2 +i3  ].patch = i1;
              pF[i+ fz+2*nLevel* i2 +i3  ].father = (int)(nPatches*(nLevel*nLevel-1)/3+i1*nLevel*nLevel+floor(i2/2)*nLevel+floor(i3/2));
              pF[i+ fz+2*nLevel* i2 +i3  ].number = pF[pF[i+ fz+2*nLevel* i2   +i3  ].father].number*4+i2%2+(i3+1)%2;
              pF[i+ fz+2*nLevel* i2 +i3  ].son[0] = nPatches*(16*nLevel*nLevel-1)/3+16*i1*nLevel*nLevel+2*i2*4*nLevel+2*i3;
              pF[i+ fz+2*nLevel* i2 +i3  ].son[1] = nPatches*(16*nLevel*nLevel-1)/3+16*i1*nLevel*nLevel+2*i2*4*nLevel+2*i3+1;
              pF[i+ fz+2*nLevel* i2 +i3  ].son[2] = nPatches*(16*nLevel*nLevel-1)/3+16*i1*nLevel*nLevel+4*(2*i2+1)*nLevel+2*i3+1;
              pF[i+ fz+2*nLevel* i2 +i3  ].son[3] = nPatches*(16*nLevel*nLevel-1)/3+16*i1*nLevel*nLevel+4*(2*i2+1)*nLevel+2*i3;
              pF[i+ fz+2*nLevel* i2 +i3  ].index_s = i3;
              pF[i+ fz+2*nLevel* i2 +i3  ].index_t = i2;
            }
            if ((i2 < 2*nLevel) && (0 < i3)){
              pElementList[fz+2*nLevel* i2   +i3-1][1] = nNodes;
              pF[i+fz+2*nLevel* i2   +i3-1].vertex[1] = nNodes;
            }
            if ((0 < i2)  && (0 < i3)){
              pElementList[fz + 2 * nLevel * (i2 - 1) + i3 - 1][2] = nNodes;
              pF[i+fz + 2 * nLevel * (i2 - 1) + i3 - 1].vertex[2] = nNodes;
            }
            if ((0 < i2)  && (i3 < 2*nLevel)){
              pElementList[fz+2*nLevel*(i2-1)+i3  ][3] = nNodes;
              pF[i+fz+2*nLevel*(i2-1)+i3  ].vertex[3] = nNodes;
            }
            ++(nNodes);
          }
        }
      }
      fz += 4*nLevel*nLevel;
    }
  }

  // loop unify - calculate circumscribed sphere on finest level
  i = nPatches*(n*n-1)/3;
  for (unsigned int l = 0; l < nPatches*n*n; ++l) {
    pF[i+l].son[0] = -1;
    pF[i+l].son[1] = -1;
    pF[i+l].son[2] = -1;
    pF[i+l].son[3] = -1;
    unify(&d1, &r1, nodeList[pF[i+l].vertex[0]], 0, nodeList[pF[i+l].vertex[2]], 0);
    unify(&d2, &r2, nodeList[pF[i+l].vertex[1]], 0, nodeList[pF[i+l].vertex[3]], 0);
    unify(&(pF[i+l].midpoint), &(pF[i+l].radius), d1, r1, d2, r2);
  }

  // assigning recursively balls from son elements
  --i;
  for(int index = i; index >=0; --index){
    unify(&d1, &r1, pF[pF[index].son[0]].midpoint, pF[pF[index].son[0]].radius,
        pF[pF[index].son[2]].midpoint, pF[pF[index].son[2]].radius);
    unify(&d2, &r2, pF[pF[index].son[1]].midpoint, pF[pF[index].son[1]].radius,
        pF[pF[index].son[3]].midpoint, pF[pF[index].son[3]].radius);
    unify(&(pF[index].midpoint), &(pF[index].radius), d1, r1, d2, r2);
  }

  return(nNodes);
}

/**
 * @brief{calculates the compression pattern. which integrals need to be
 * computed and which can be left out of the system matrix calculation}
 */
unsigned int GenericAnsatzFunction :: compression(SparseMatrix *T){
  double maxRadius = 0.0;
  double **c1, **c2;
  double d1, d2;
  unsigned int m1, m2;
  unsigned int ind1, ind2;
  unsigned int nnz;
  unsigned int *rn;
  FILE* debugFile;

  for(unsigned int i = 0; elementTree.element[i].level == 0; ++i)
    if( maxRadius < elementTree.element[i].radius)
      maxRadius = elementTree.element[i].radius;
  c1 = (double**) calloc(1,sizeof(double*)*(nLevels+1));
  c2 = (double**) calloc(1,sizeof(double*)*(nLevels+1));
  for(unsigned int i = 0; i <= nLevels; ++i){
    c1[i] = (double*) calloc(1,sizeof(double)*(i+1));
    c2[i] = (double*) calloc(1,sizeof(double)*(i+1));
    for(unsigned int j = 0; j <= i; ++j){
      c1[i][j] = a*pow(2,(nLevels*(2*dp-op)-(i+j)*(dp+td))/(2*td+op));
      d1 = pow(2,(nLevels*(2*dp-op)-(i+j)*dp-i*td)/(td+op));    // alter Parameter
      d2 = pow(2,(2*nLevels-(i+j))*(2*dp-op)/((2*td+op)*(td+op)))*pow(2,(nLevels*(2*dp-op)-i*(dp+td+1)-j*(dp-1))/(td+op));
      if (d1 > d2) c2[i][j] = a*d1; else c2[i][j] = a*d2;
      if (c1[i][j] < a/(1<<j)) c1[i][j] = a/(1<<j);
      if (c2[i][j] < a/(1<<i)) c2[i][j] = a/(1<<i);
      if (c1[i][j] > c2[i][j]) c1[i][j] = c2[i][j];
      if ((i < minLevel) || (j < minLevel)) c1[i][j] = c2[i][j];
      c1[i][j] *= maxRadius/scalingFactor;              // Gebiet relativieren
    }
  }
#ifdef DEBUG
  debugFile = fopen("debug.out","a");
  fprintf(debugFile,">>> COMP CC %lf %lf %lf\n", td, dp, op);
  for(unsigned int i  = 0; i < nLevels+1; ++i){
    for(unsigned int j = 0; j < i+1; ++j){
      fprintf(debugFile, "[%lf %lf]", c1[i][j], c2[i][j]);
    }
    fprintf(debugFile,"\n");
  }
  fprintf(debugFile,"\n<<< COMP CC\n");
  fflush(debugFile);
  fclose(debugFile);
#endif
  initPattern2(T,waveletList.sizeWaveletList,waveletList.sizeWaveletList,20);

  computeBoundingBoxes();
#ifdef DEBUG
  debugFile = fopen("debug.out", "a");
  fprintf(debugFile,">>> COMPRESSION BB\n"); 
  fprintf(debugFile, "%d\n", waveletList.sizeWaveletList);
  for(unsigned int i  = 0; i < waveletList.sizeWaveletList; ++i){
    fprintf(debugFile, "%lf %lf %lf - %lf %lf %lf\n",B[i].rx, B[i].ry, B[i].rz, B[i].mx, B[i].my, B[i].mz);
  }
  fprintf(debugFile,"<<< COMPRESSION BB\n");
  fclose(debugFile);
#endif

  for(unsigned int i = 0; i<waveletList.sizeWaveletList && waveletList.W[i].level < minLevel; ++i){
    for(unsigned int j = 0; j <=i; ++j)
      setPattern(T,i,j);
  }
#ifdef DEBUG
  debugFile = fopen("debug.out", "a");
  fprintf(debugFile,">>> COMPRESSION CC");  
  for(unsigned int i  = 0; i < T->m; ++i){
    fprintf(debugFile,"\n%d %d\n", T->row_number[i], minLevel);
    for(unsigned int j = 0; j < T->row_number[i]; ++j){
      fprintf(debugFile, "%d ", T->index[i][j]);
    }
  }
  fprintf(debugFile,"\n<<< COMPRESSION CC\n");
  fflush(debugFile);
  fclose(debugFile);
#endif

#ifdef DEBUG
  debugFile = fopen("debug.out", "a");
  fprintf(debugFile,">>> WAVELETWAVELETCRITERION\n"); 
  fclose(debugFile);
#endif
#ifdef DEBUG
  debugFile = fopen("debug.out", "a");
#endif
  for(unsigned int i = 0; waveletList.W[i].level < nLevels; ++i){
    m1 = waveletList.W[i].level+1;
  // bestimme aus den Bloecken (m1-1,1:m1-1) die Bloecke (m1,1:m1)
    for(unsigned int j = 0; j < T->row_number[i]; ++j){
      ind2 = T->index[i][j];
      m2 = waveletList.W[ind2].level;
      if(m1 == m2+1){
        for(unsigned int k = 0; k < waveletList.W[i].noSons; ++k){
          ind1 = waveletList.W[i].son[k];
          //if ((ind2 <= ind1) && (waveletWaveletCriterion(ind1,ind2,c1[m1][m2],c2[m1][m2]))){
          if(waveletWaveletCriterion(ind1, ind2, c1[m1][m2], c2[m1][m2])){ 
            setPattern(T, ind1, ind2);
#ifdef DEBUG
  fprintf(debugFile,"set 1 %d %d\n", ind1, ind2);
#endif
            for(unsigned int l = 0; l < waveletList.W[ind2].noSons; ++l){
#ifdef DEBUG
  fprintf(debugFile,"set 2 %d %d\n",ind1, waveletList.W[ind2].son[l]);
#endif
              if ((waveletList.W[ind2].son[l] <= ind1) && (waveletWaveletCriterion(ind1,waveletList.W[ind2].son[l],c1[m1][m1],c2[m1][m1]))){  
              //if(waveletWaveletCriterion(ind1, waveletList.W[ind2].son[l], c1[m1][m1], c2[m1][m1])) {
                setPattern(T, ind1, waveletList.W[ind2].son[l]);
#ifdef DEBUG
  fprintf(debugFile,"set %d %d\n",ind1, waveletList.W[ind2].son[l]);
#endif
              }
            }
          }
        }
      // wegen halbem Diagonalblock
        for(unsigned int k = 0; k < waveletList.W[ind2].noSons;++k){
          if(waveletWaveletCriterion(waveletList.W[ind2].son[k], i, c1[m1][m2], c2[m1][m2])){
            setPattern(T, waveletList.W[ind2].son[k], i);
#ifdef DEBUG
  fprintf(debugFile,"set 3 %d %d\n", waveletList.W[ind2].son[k], i);
#endif
            for(unsigned int l = 0; l < waveletList.W[i].noSons;++l){
#ifdef DEBUG
  fprintf(debugFile,"set 4 %d %d\n", waveletList.W[ind2].son[k], waveletList.W[i].son[l]);
#endif
              if ((waveletList.W[i].son[l] <= waveletList.W[ind2].son[k]) && (waveletWaveletCriterion(waveletList.W[ind2].son[k],waveletList.W[i].son[l],c1[m1][m1],c2[m1][m1]))){  
              //if (waveletWaveletCriterion(waveletList.W[ind2].son[k], waveletList.W[i].son[l], c1[m1][m1], c2[m1][m1])){
                setPattern(T, waveletList.W[ind2].son[k], waveletList.W[i].son[l]);
#ifdef DEBUG
  fprintf(debugFile,"set %d %d\n", waveletList.W[ind2].son[k], waveletList.W[i].son[l]);
#endif
              }
            }
          }
        }
      } else {
        for(unsigned int k = 0; k < waveletList.W[i].noSons; ++k){
          ind1 = waveletList.W[i].son[k];
#ifdef DEBUG
  fprintf(debugFile,"set 5 %d %d\n",ind1, ind2);
#endif
          //if ((ind2 <= ind1) && (waveletWaveletCriterion(ind1,ind2,c1[m1][m2],c2[m1][m2]))){  
          if(waveletWaveletCriterion(ind1, ind2, c1[m1][m2], c2[m1][m2])){
            setPattern(T, ind1, ind2);
#ifdef DEBUG
  fprintf(debugFile,"set %d %d\n",ind1, ind2);
#endif
          }
        }
      }
    }
  }
#ifdef DEBUG
  fclose(debugFile);
#endif
#ifdef DEBUG
  debugFile = fopen("debug.out", "a");
  fprintf(debugFile,"<<< WAVELETWAVELETCRITERION\n"); 
  fclose(debugFile);
#endif

  nnz = 0;
  rn = (unsigned int*) calloc(waveletList.sizeWaveletList, sizeof(unsigned int));
  for(unsigned int i = 0; i < waveletList.sizeWaveletList; ++i){
    rn[i] += T->row_number[i];
    for(unsigned int j = 0; j+1 < T->row_number[i]; ++j){
      rn[T->index[i][j]]++;
    }
  }
  for(unsigned int i = 0; i < waveletList.sizeWaveletList; ++i){
    T->index[i] = (unsigned int*)realloc(T->index[i], (rn[i]+1)*sizeof(unsigned int));
    T->max_row_number[i] = rn[i];
    T->index[i][rn[i]] = waveletList.sizeWaveletList;
    nnz = nnz+rn[i];
  }
  for(unsigned int i = 0; i < waveletList.sizeWaveletList; ++i){
    for(unsigned int j = 0; j+1 < T->row_number[i]; ++j){
      T->index[T->index[i][j]][T->row_number[T->index[i][j]]++] = i;
    }
  }
#ifdef DEBUG
  debugFile = fopen("debug.out", "a");
  fprintf(debugFile,">>> COMPRESSION C");  
  for(unsigned int i  = 0; i < T->m; ++i){
    fprintf(debugFile,"\n%d\n", T->row_number[i]);
    for(unsigned int j = 0; j < T->row_number[i]; ++j){
      fprintf(debugFile, "%d ", T->index[i][j]);
    }
  }
  fprintf(debugFile,"\n<<< COMPRESSION C\n");
  fclose(debugFile);
#endif
  printf("A-priori compression:            %.5f %%\n",100.0*nnz/waveletList.sizeWaveletList/waveletList.sizeWaveletList);
  finishPattern2(T);
#ifdef DEBUG
  debugFile = fopen("debug.out", "a");
  fprintf(debugFile,">>> FINISH");  
  for(unsigned int i  = 0; i < T->m; ++i){
    fprintf(debugFile,"\n%d\n", T->row_number[i]);
    for(unsigned int j = 0; j < T->row_number[i]; ++j){
      fprintf(debugFile, "%d ", T->index[i][j]);
    }
  }
  fprintf(debugFile,"\n<<< FINISH\n");
  fclose(debugFile);
#endif
  free (rn);
  for(unsigned int i = 0 ; i <= nLevels; ++i){
    free(c1[i]);
    free(c2[i]);
  }
  free (c1);
  free (c2);
  freeBoundingBoxes();
  return 0;
}

unsigned int GenericAnsatzFunction :: postProc(SparseMatrix *T){
  double		**c;        // cut-off parameter
  double		*D1, *D2;   // diagonal entries of the systemmatrix
  double		*v1, *v2;   // vector for new values
  unsigned int	*w;     // vector for new indeces
  unsigned int    nnz;  // number of non-zero elements of systemmatrix
  unsigned int    k;

  // calculate diagonal entries for value1 and value2
  D1 = (double*) malloc(T->n*sizeof(double));
  D2 = (double*) malloc(T->n*sizeof(double));
  for (unsigned int i=0; i<T->n; ++i) {
    getSparse2(T,i,i,&D1[i],&D2[i]);
    D1[i] = sqrt(fabs(D1[i]));
    D2[i] = sqrt(fabs(D2[i]));
  }

  // calculate cut-off parameter
  c = (double**) malloc((nLevels+1)*sizeof(double*));
  for (int i=0; i<= (int)nLevels; ++i){
    c[i] = (double*) malloc((nLevels+1)*sizeof(double));
    for (int j=0; j<= (int)nLevels; ++j){
      // LinAF 
      c[i][j] = pow(0.5,(nLevels-0.5*(i+j))*(2*dp-op)/(2*td+op)); // false?!
      //c[i][j] = pow(0.5,(2*nLevels-(i+j))*(2*dp-op)/(2*td+op));
      if (c[i][j] > pow(0.5,fabs(i-j))) c[i][j] = pow(0.5,fabs(i-j));
      //c[i][j] *= b * pow(0.5,(dp-0.5*op)*(2*nLevels-(i+j)));// added -op*nLevels
      c[i][j] *= b * pow(0.5,(dp)*(2*nLevels-(i+j)));
    }
  }

  // a-posteriori compression
  nnz = 0;
  for (unsigned int i=0; i<T->n; ++i){  
    // find relevant matrix entries
    k=0;
    v1 = (double*) malloc(T->row_number[i]*sizeof(double));
    v2 = (double*) malloc(T->row_number[i]*sizeof(double));
    w = (unsigned int*) malloc((T->row_number[i]+1)*sizeof(unsigned int));

    for (unsigned int j=0; j<T->row_number[i]; ++j){
      if( (fabs(T->value1[i][j]) >= D1[i]*D1[T->index[i][j]]*c[waveletList.W[i].level][waveletList.W[T->index[i][j]].level]) ||
          (fabs(T->value2[i][j]) >= D2[i]*D2[T->index[i][j]]*c[waveletList.W[i].level][waveletList.W[T->index[i][j]].level]) ){
        v1[k] = T->value1[i][j];
        v2[k] = T->value2[i][j];
        w[k]  = T->index[i][j];
        ++k;
      }
    }

    // copy new values and free memory for old ones
    free(T->index[i]);
    free(T->value1[i]);
    free(T->value2[i]);
    T->max_row_number[i] = T->row_number[i] = k;
    T->value1[i] = (double*)       realloc(v1, k   *sizeof(double));
    T->value2[i] = (double*)       realloc(v2, k   *sizeof(double));
    T->index[i] = (unsigned int*) realloc(w,(k+1)*sizeof(unsigned int));
    T->index[i][k] = T->n;
    nnz += k;
  }
   
  // free memory of function
  //printf("A-posteriori compression:        %.5f %%\n",100.0*nnz/T->n/T->n);
  for (unsigned int i=0; i<=nLevels; ++i) free(c[i]);
  free(D1);
  free(D2);
  free(c);
  return 100.0*nnz/T->n/T->n;
}

void GenericAnsatzFunction :: elementElementInteraction(double *c, unsigned int ind1, unsigned int ind2, double prec, Cubature * Q, double SingleLayer(Vector3, Vector3), double DoubleLayer(Vector3, Vector3, Vector3), double Identity){
  signed int j;
  double dist;
  signed int g1, g2;
  unsigned int s, t;
  unsigned int CASE;
  double *a;
  unsigned int sizeC = 3*noPhi*noPhi;
  
  a = (double*) calloc(sizeC*4,sizeof(double));
  memset(c, 0, sizeC*sizeof(double));
  
  if (ind1 >= ind2){
    j = searchIntegral(&elementTree.element[ind1].interaction, ind2);
    if(j!=-1){
      //integral already computed, copy value
      memcpy(c, elementTree.element[ind1].interaction.value[j], sizeC*sizeof(double));
      free(a);
      return;
    }
  } else {
    j = searchIntegral(&elementTree.element[ind2].interaction, ind1);
    if (j != -1) {
      permutate(c, elementTree.element[ind2].interaction.value[j]);
      free(a);
      return;
    }
  }

  // calculate distance between two elements
  dist = distance(ind1, ind2);
  // quadrature with precision prec
  //printf("%lf %d %d\n", dist, ind1, ind2);
  if (elementTree.element[ind1].level == elementTree.element[ind2].level) {
    // both elements are on the same level no need for division
    // nothing in common
    if(dist > eps) CASE = 1; 
    // same patch
    else if( ind1 == ind2 ) CASE = 2; 
    // common vertex or edge
    else CASE = compare (ind1, ind2, &s, &t);

    // quadrature with precision prec
    if(dist *(1<<elementTree.element[ind1].level) < 1) dist = 1./(1<<elementTree.element[ind1].level);
    quadratureGrade(&g1, &g2, elementTree.element[ind1].level, elementTree.element[ind2].level, dist, prec); // TODO

    // choose the integration routine
    //printf("%d %d %d %lf %d %d\n", CASE, g1, g2, dist, ind1, ind2);
    switch(CASE) {
      case 1:
        integrateNoProblem(c, ind1, ind2, g1, &Q[g1], &Q[g1], SingleLayer, DoubleLayer);
        //memset(c, 1, sizeof(double)*sizeC);
        break;
      case 2:
        integratePatch(c, ind1, &Q[g1], SingleLayer, DoubleLayer, Identity);
        //memset(c, 2, sizeof(double)*sizeC);
        break;
      case 3:
        integrateEdge(c, ind1, ind2, s, t, &Q[g1], SingleLayer, DoubleLayer, Identity);
        //memset(c, 3, sizeof(double)*sizeC);
        break;
      case 4:
        integratePoint(c, ind1, ind2, s, t, &Q[g1], SingleLayer, DoubleLayer, Identity);
        //memset(c, 4, sizeof(double)*sizeC);
        break;
      default:
        //memset(c,0,sizeC*sizeof(double));
        printf("ERROR : CASE!=1,2,3,4\n");
    }
  //printf("%lf %lf %lf\n",c[0], c[1], c[2]);
  } else {
    if(dist*(1<<elementTree.element[ind2].level) >=q ){
      quadratureGrade(&g1, &g2, elementTree.element[ind1].level, elementTree.element[ind2].level, dist, prec); // TODO
      integrateNoProblem(c, ind1, ind2, g1, &Q[g1], &Q[g2], SingleLayer, DoubleLayer);
      //memset(c, 1, sizeof(double)*sizeC);
      //printf("%d %d %d %d\n", ind1, ind2, g1, g2);
      //for(unsigned int k = 0; k < 48; ++k) printf("%lf ", c[k]);
      //printf("\n");
    } else {
      //divide element
      elementElementInteraction(&a[0*sizeC], ind1, elementTree.element[ind2].son[0], prec, Q, SingleLayer, DoubleLayer, Identity);
      elementElementInteraction(&a[1*sizeC], ind1, elementTree.element[ind2].son[1], prec, Q, SingleLayer, DoubleLayer, Identity);
      elementElementInteraction(&a[2*sizeC], ind1, elementTree.element[ind2].son[2], prec, Q, SingleLayer, DoubleLayer, Identity);
      elementElementInteraction(&a[3*sizeC], ind1, elementTree.element[ind2].son[3], prec, Q, SingleLayer, DoubleLayer, Identity);
      //for(int lr = 0 ; lr < sizeC*4; ++lr) printf("%lf\n", a[lr]);
      calculateIntegral(a, c);
      //memset(c, 1, sizeof(double)*sizeC);
    }
  }
  //printf("%lf %lf %lf\n",c[0], c[1], c[2]);

  // save integrals for further use
  if (ind1 >= ind2) setIntegral(&elementTree.element[ind1].interaction, ind2, c);
  else {
    permutate(a,c);
    setIntegral(&elementTree.element[ind2].interaction, ind1, a);
  }
  free(a);
  return;
}

signed int GenericAnsatzFunction::searchIntegral(intvector *I, unsigned int i){
  unsigned int low, mid, high;

  if (I->integralNumber == 0) return -1;

  // binary search
  low = 0;
  high = I->integralNumber-1;

  while(low < high){
    mid = (low + high)/2;
    if( I->index[mid] < i) low = mid+1;
    else if (I->index[mid] > i) high = mid;
    else return mid;
  }

  if (I->index[low] == i) return low;
  else return -1;
}

void GenericAnsatzFunction::setIntegral(intvector *I, unsigned int i, double* z){
  // der Eintrag darf nicht vorhanden sein
  unsigned int rn, low, mid, high;

  unsigned int sizeC = 3*noPhi*noPhi;
  rn = I->integralNumber++;
  if(rn%delta == 0){
    if(rn == 0){
      I->index = (unsigned int *) malloc(delta*sizeof(unsigned int));
      I->value = (double **) malloc(delta*sizeof(double*));
      I->value[0] = (double *) malloc(sizeC*sizeof(double));
      I->index[0] = i;
      memcpy(I->value[0],z,sizeC*sizeof(double));
      return;
    } else {
      I->index = (unsigned int*) realloc(I->index,(rn+delta)*sizeof(unsigned int));
      I->value = (double**) realloc(I->value, (rn+delta)*sizeof(double*));
    }
  }

  // can we simply attach integral value
  if(I->index[rn-1]< i){ 
    I->index[rn] = i;
    I->value[rn] = (double*) malloc(sizeC*sizeof(double));
    memcpy(I->value[rn], z, sizeC*sizeof(double));
    return;
  }

  // search entry
  //I->value[rn] = (double*) malloc(sizeC*sizeof(double));
  low = 0;
  high = rn-1;
  while (low< high){
    mid = (low + high)/2;
    if(I->index[mid] < i) low = mid+1;
    else high = mid;
  }

  //add to list
  memmove(&I->index[low+1], &I->index[low], (rn-low)*sizeof(unsigned int));
  memmove(&I->value[low+1], &I->value[low], (rn-low)*sizeof(double*));
  I->index[low] = i;
  I->value[low] = (double*) malloc(sizeC*sizeof(double));
  memcpy(I->value[low],z,sizeC*sizeof(double));
  return;
}

void GenericAnsatzFunction :: inv_A_times_x(double *x){
  double u, uStart, v, omg;
  double *r, *d, *Ad, *z;

  r  = (double*)calloc(waveletList.sizeWaveletList,sizeof(double));
  d  = (double*)malloc(sizeof(double)*waveletList.sizeWaveletList);
  Ad = (double*)malloc(sizeof(double)*waveletList.sizeWaveletList);
  z  = (double*)malloc(sizeof(double)*waveletList.sizeWaveletList);

  // r = x - T*A*T'*x, d = r und u = (r,r)
  u = 0;
  memcpy(d,x,waveletList.sizeWaveletList*sizeof(double));
  tdwt(d);
  for (unsigned int i=0; i<G->n; ++i){
    for (unsigned int j=0; j<G->row_number[i]; ++j) r[i] -= G->value[i][j] * d[G->index[i][j]];
  }
  dwt(r);
  for (unsigned int i=0; i<waveletList.sizeWaveletList; ++i){
    r[i] += x[i];
    d[i] = r[i];
    u += r[i] * r[i];
  }
  uStart = u;

  // Iteration
  while (sqrt(u/uStart)>eps){
    // Ad = A*d
    memcpy(z,d,waveletList.sizeWaveletList*sizeof(double));
    memset(Ad,0,waveletList.sizeWaveletList*sizeof(double));
    tdwt(z);
    for (unsigned int i=0; i<G->n; ++i){
      for (unsigned int j=0; j<G->row_number[i]; ++j) Ad[i] += G->value[i][j] * z[G->index[i][j]];
    }
    dwt(Ad);

    // omg = u / (d,Ad) und v = u
    omg = 0;
    for (unsigned int i=0; i<waveletList.sizeWaveletList; ++i) omg += d[i] * Ad[i];
    omg = u/omg;
    v = u;

    // x = x + omg * d, r = r - omg * Ad und u = (r,r)
    u = 0;
    for (unsigned int i=0; i<waveletList.sizeWaveletList; ++i){
      x[i] += omg * d[i];
      r[i] -= omg * Ad[i];
      u += r[i] * r[i];
    }

    // d = r +  u/v * d
    omg = u/v;
    for (unsigned int i=0; i<waveletList.sizeWaveletList; ++i) d[i] = r[i] + omg * d[i];
  }

  free(Ad);
  free(z);
  free(r);
  free(d);
  return;
}

void GenericAnsatzFunction ::  precond(double *a, double *b){
  // apply Wavelet preconditioner
  for(unsigned int i = 0; i < waveletList.sizeWaveletList; ++i) a[i] = pow(2, -0.5*op*waveletList.W[i].level)*b[i];

  inv_A_times_x(a);

  // apply Wavelet preconditioner
  for(unsigned int i = 0; i < waveletList.sizeWaveletList; ++i) a[i] *= pow(2, -0.5*op*waveletList.W[i].level);

}

double GenericAnsatzFunction::distance(int e1, int e2){
  double dx, dy, dz;
  double c1, c2;
  double dist;

  dx = elementTree.element[e1].midpoint.x - elementTree.element[e2].midpoint.x;
  dy = elementTree.element[e1].midpoint.y - elementTree.element[e2].midpoint.y;
  dz = elementTree.element[e1].midpoint.z - elementTree.element[e2].midpoint.z;

  dist = sqrt(dx*dx + dy*dy + dz*dz) - elementTree.element[e1].radius - elementTree.element[e2].radius;

  c1 = elementTree.element[e1].radius*(1<<elementTree.element[e1].level);
  c2 = elementTree.element[e2].radius*(1<<elementTree.element[e2].level);

  if(c1 < c2) return scalingFactor*dist/c2;
  else return scalingFactor*dist/c1;
}

int GenericAnsatzFunction::print_geometry(double* rho, char* dname) {
  unsigned int i = 0;
  unsigned int j = 0;
  unsigned int max = 0;
  FILE *f = NULL;
  Vector3 auxPoint;
  int n = 2<<nLevels;
  int index = nPatches*(n*n/4-1)/3;
  et_node *Ef = elementTree.element;

  // get length of point list P
  max = nNodes;

  // visualize the surface of the geometry in vtk format
  f = fopen(dname, "w");

  // vtk format header
  fprintf(f, "# vtk DataFile Version 3.1\n");
  fprintf(f, "this file hopefully represents my surface now\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  // print point list
  fprintf(f, "POINTS %d FLOAT\n", max);
  for (i = 0; i < max; ++i) {
    fprintf(f, "%20.16f\t%20.16f\t%20.16f\n", nodeList[i].x,
        nodeList[i].y, nodeList[i].z);
  }
  fprintf(f, "\n");

  // print element list
  fprintf(f, "CELLS %d %d\n", nFunctions, 5 * nFunctions);
  for (i = 0; i < nFunctions; ++i)
    fprintf(f, "%d %d %d %d %d\n", 4, Ef[index+i].vertex[0], Ef[index+i].vertex[1],
        Ef[index+i].vertex[2], Ef[index+i].vertex[3]);
  fprintf(f, "\n");

  fprintf(f, "CELL_TYPES %d\n", nFunctions);
  for (i = 0; i < nFunctions; ++i)
    fprintf(f, "%d\n", 9);
  fprintf(f, "\n");

  // print z-values of the geometry and solved density for visualization
  fprintf(f, "POINT_DATA %d\n", max);
  fprintf(f, "SCALARS Z-value FLOAT\n");
  fprintf(f, "LOOKUP_TABLE default\n");
  if(rho != NULL){
    for (i = 0; i < max; ++i)
      fprintf(f, "%20.16f\n", rho[i]);
  } else {
    for (i = 0; i < max; ++i)
      fprintf(f, "%20.16f\n", nodeList[i].z);
  }
  fprintf(f, "\n");

  fprintf(f, "CELL_DATA %d\n", nFunctions);
  fprintf(f, "SCALARS Patch FLOAT\n");
  fprintf(f, "LOOKUP_TABLE default\n");
  for (i = 0; i < nFunctions; ++i)
    fprintf(f, "%d\n", Ef[index+i].patch);
  fprintf(f, "\n");
  
  Vector2 xi;
  Vector3 normal;
  double ht = 1./(1<<Ef[index].level);
  fprintf(f, "NORMALS Cell_Norm FLOAT\n");
  fprintf(f, "LOOKUP_TABLE default\n");
  if (rho != NULL) {
    for (i = 0; i < nFunctions; ++i)
      fprintf(f, "%20.16f\n", rho[i]);
  }else{
    for(i = 0; i < nFunctions; ++i){
      xi.x = ht*Ef[index+i].index_s;
      xi.y = ht*Ef[index+i].index_t;
      normal = interCoeff->n_Chi(xi,Ef[index+i].patch);
      normal = vector3SMul(1./vector3Norm(normal),normal);
      fprintf(f, "%20.16f %20.16f %20.16f\n",normal.x, normal.y, normal.z);
    }
  }

  fclose(f);
  return 0;
}

int GenericAnsatzFunction::printElement(char* dname, int element) {
#ifdef PRINT_DEBUG
  // print out the first element!
  int i = 0;
  int j = 0;
  int i1, i2, i3, j1, j2, j3;
  int max = 0;
  FILE *f = NULL;
  int pointsToAdd = 1;
  Vector3 auxPoint;
  int n = 1<<nLevels;
  int index = nPatches*(n*n-1)/3;
  et_node *Ef = elementTree.element;

  max = (2 + pointsToAdd) * (2 + pointsToAdd);
  double w[2 + pointsToAdd];
  for (i = 1; i <= 2 + pointsToAdd; ++i)
    w[i] = i / (2 + pointsToAdd);

  // visualize the surface of the geometry in vtk format              
  f = fopen(dname, "w");

  // vtk format header                              
  fprintf(f, "# vtk DataFile Version 3.1\n");
  fprintf(f, "this file hopefully represents my surface now\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  // print point list                               
  fprintf(f, "POINTS %d FLOAT\n", max);

  //use the fact that i know that i add only one point :)
  double h = 1. / (1 << nLevels);
  i1 = Ef[index+element].patch;
  i2 = Ef[index+element].index_t;
  i3 = Ef[index+element].index_s;
  i = 0;
  fprintf(f, "%20.16f\t%20.16f\t%20.16f\n",
      (nodeList[Ef[index+element].vertex[i]]).x,
      (nodeList[Ef[index+element].vertex[i]]).y,
      (nodeList[Ef[index+element].vertex[i]]).z);
  //auxPoint = interCoeff->Chi( Vector2(h*Ef[element].index_s, h*Ef[element].index_t), Ef[element].patch); 
  //fprintf(f, "%20.16f\t%20.16f\t%20.16f\n", auxPoint.x, auxPoint.y, auxPoint.z);

  //printf("%f %f %f %f %f %f\n", h*Ef[element].index_s, h*Ef[element].index_s+h*0.5, h*Ef[element].index_s+h, h*Ef[element].index_t, h*Ef[element].index_t+0.5*h, h*Ef[element].index_t+h);
  //printf("%f %f %f %f %f %f\n", h*Ef[0].index_s, h*Ef[0].index_s+h*0.5, h*Ef[0].index_s+h, h*Ef[0].index_t, h*Ef[0].index_t+0.5*h, h*Ef[0].index_t+h);
  auxPoint = interCoeff->Chi(
      Vector2(h * Ef[index+element].index_s + 0.5 * h,
          h * Ef[index+element].index_t), Ef[index+element].patch);
  fprintf(f, "%20.16f\t%20.16f\t%20.16f\n", auxPoint.x, auxPoint.y,
      auxPoint.z);

  i = 1;
  fprintf(f, "%20.16f\t%20.16f\t%20.16f\n",
      (nodeList[Ef[index+element].vertex[i]]).x,
      (nodeList[Ef[index+element].vertex[i]]).y,
      (nodeList[Ef[index+element].vertex[i]]).z);
  //auxPoint = interCoeff->Chi( Vector2(h*Ef[element].index_s+h, h*Ef[element].index_t), Ef[element].patch); 
  //fprintf(f, "%20.16f\t%20.16f\t%20.16f\n", auxPoint.x, auxPoint.y, auxPoint.z);

  auxPoint = interCoeff->Chi(
      Vector2(h * Ef[index+element].index_s,
          h * Ef[index+element].index_t + 0.5 * h), Ef[index+element].patch);
  fprintf(f, "%20.16f\t%20.16f\t%20.16f\n", auxPoint.x, auxPoint.y,
      auxPoint.z);

  auxPoint = interCoeff->Chi(
      Vector2(h * Ef[index+element].index_s + 0.5 * h,
          h * Ef[index+element].index_t + 0.5 * h), Ef[index+element].patch);
  fprintf(f, "%20.16f\t%20.16f\t%20.16f\n", auxPoint.x, auxPoint.y,
      auxPoint.z);

  auxPoint = interCoeff->Chi(
      Vector2(h * Ef[index+element].index_s + h,
          h * Ef[index+element].index_t + 0.5 * h), Ef[index+element].patch);
  fprintf(f, "%20.16f\t%20.16f\t%20.16f\n", auxPoint.x, auxPoint.y,
      auxPoint.z);
  i = 3;
  fprintf(f, "%20.16f\t%20.16f\t%20.16f\n",
      (nodeList[Ef[index+element].vertex[i]]).x,
      (nodeList[Ef[index+element].vertex[i]]).y,
      (nodeList[Ef[index+element].vertex[i]]).z);
  //auxPoint = interCoeff->Chi( Vector2(h*Ef[element].index_s, h*Ef[element].index_t+h), Ef[element].patch); 
  //fprintf(f, "%20.16f\t%20.16f\t%20.16f\n", auxPoint.x, auxPoint.y, auxPoint.z);

  auxPoint = interCoeff->Chi(
      Vector2(h * Ef[element].index_s + 0.5 * h,
          h * Ef[element].index_t + h), Ef[index+element].patch);
  fprintf(f, "%20.16f\t%20.16f\t%20.16f\n", auxPoint.x, auxPoint.y,
      auxPoint.z);

  i = 2;
  fprintf(f, "%20.16f\t%20.16f\t%20.16f\n",
      (nodeList[Ef[index+element].vertex[i]]).x,
      (nodeList[Ef[index+element].vertex[i]]).y,
      (nodeList[Ef[index+element].vertex[i]]).z);
  //auxPoint = interCoeff->Chi( Vector2(h*Ef[element].index_s+h, h*Ef[element].index_t+h), Ef[element].patch); 
  //fprintf(f, "%20.16f\t%20.16f\t%20.16f\n", auxPoint.x, auxPoint.y, auxPoint.z);
  fprintf(f, "\n");

  // print element list                             
  fprintf(f, "CELLS %d %d\n", 1, 5 * 1);
  fprintf(f, "%d %d %d %d %d\n", 4, Ef[index+element].vertex[0],
      Ef[index+element].vertex[1], Ef[index+element].vertex[2],
      Ef[index+element].vertex[3]);
  fprintf(f, "\n");

  fprintf(f, "CELL_TYPES %d\n", 1);
  fprintf(f, "%d\n", 9);
  fprintf(f, "\n");

  // print z-values of the geometry and solved density for visualization      
  fprintf(f, "POINT_DATA %d\n", max);
  fprintf(f, "SCALARS Z-value FLOAT\n");
  fprintf(f, "LOOKUP_TABLE default\n");
  for (i = 0; i < 4; ++i)
    fprintf(f, "%20.16f\n", (nodeList[Ef[index+element].vertex[i]]).z);
  fprintf(f, "\n");

  fprintf(f, "CELL_DATA %d\n", 1);
  fprintf(f, "SCALARS Cell_Density FLOAT\n");
  fprintf(f, "LOOKUP_TABLE default\n");

  fclose(f);
#endif
  return 0;
}

/**
 * @brief { checks two elements with index e1 and e2 for a common edge (return
 * 3) or vertex (return 4), ind1 and ind2 are the corresponding rotations of
 * the elements. In case of no common intersection returns 1. In case of same
 * patch, the comparison is done via indeces. In case of different patches the
 * comparison is done via differences}
 *
 * @param e1,e2 index of elements to compare
 * @param ind1, ind2 indeces for the rotations
 */
unsigned int GenericAnsatzFunction::compare(unsigned int e1, unsigned int e2,
    unsigned int *ind1, unsigned int *ind2) {
  unsigned int *F1 = elementTree.element[e1].vertex;
  unsigned int *F2 = elementTree.element[e2].vertex;
  if(elementTree.element[e1].patch == elementTree.element[e2].patch){
    // same patch
    for (*ind1 = 0; *ind1 < 4; ++(*ind1)){
      //common point?
      for( *ind2 = 0; *ind2 < 4; ++(*ind2)){
        if( F1[*ind1] == F2[*ind2] ){
          //common edge?
          //special case
          if( F1[3] == F2[(*ind2+1)%4]){
            *ind1 = 3;
            return 3;
          } else if (F1[(*ind1+1)%4] == F2[(*ind2+3)%4]) {
					  // normal case
            *ind2 = (*ind2+3)%4;
            return(3);
          }
          else return(4);
        }
      }
    }
    return 1;
  }else{
    //different patches, calculate difference vector
    Vector3 d;

    // search common point
    for (*ind1 = 0; *ind1 < 4; ++(*ind1)) {
      for (*ind2 = 0; *ind2 < 4; ++(*ind2)) {
        d.x = nodeList[F1[*ind1]].x - nodeList[F2[*ind2]].x;
        d.y = nodeList[F1[*ind1]].y - nodeList[F2[*ind2]].y;
        d.z = nodeList[F1[*ind1]].z - nodeList[F2[*ind2]].z;
        if (d.x * d.x + d.y * d.y + d.z * d.z < 1e-8) {
          // found common point, search common edge
          d.x = nodeList[F1[(*ind1 + 1) % 4]].x
            - nodeList[F2[(*ind2 + 3) % 4]].x;
          d.y = nodeList[F1[(*ind1 + 1) % 4]].y
            - nodeList[F2[(*ind2 + 3) % 4]].y;
          d.z = nodeList[F1[(*ind1 + 1) % 4]].z
            - nodeList[F2[(*ind2 + 3) % 4]].z;
          if (d.x * d.x + d.y * d.y + d.z * d.z < 1e-8) {
            *ind2 = (*ind2 + 3) % 4;// normal case: second point at ind1+1, ind1-1
            return (3);
          } else if (*ind1 == 0) {
            d.x = nodeList[F1[3]].x
              - nodeList[F2[(*ind2 + 1) % 4]].x; // special case, 1st point at ind1=0, second point can be ind2=3
            d.y = nodeList[F1[3]].y
              - nodeList[F2[(*ind2 + 1) % 4]].y; 
            d.z = nodeList[F1[3]].z
              - nodeList[F2[(*ind2 + 1) % 4]].z;
            if (d.x * d.x + d.y * d.y + d.z * d.z < 1e-8) {
              *ind1 = 3;
              return (3);
            }
          }
          return (4);
        }
      }
    }
    return (1);
  }
}

// allocate far field values
void GenericAnsatzFunction::initRandwerte(int g_max) {
  pRandWerte = (Randwerte*) malloc(g_max * sizeof(Randwerte));
  for (int i1 = 0; i1 < g_max; ++i1) {
    pRandWerte[i1].Chi = (Vector3*) malloc((i1+1)*(i1+1) * sizeof(Vector3));
    pRandWerte[i1].n_Chi = (Vector3*) malloc((i1+1)*(i1+1) * sizeof(Vector3));
    pRandWerte[i1].det_dChi = (double*) malloc((i1+1)*(i1+1) * sizeof(double));
    pRandWerte[i1].noP = 0;
  }
  return;
}

// erases the content of far fiel values
void GenericAnsatzFunction::resetRandwerte(int g_max) {
  for (int i1 = 0; i1 < g_max; ++i1) {
    pRandWerte[i1].noP = 0;
  }
  return;
}

// frees the memory allocated for far field values
void GenericAnsatzFunction::freeRandwerte() {
  for(unsigned int i1 = 0; i1 < 10; ++i1){
    free(pRandWerte[i1].Chi);
    free(pRandWerte[i1].n_Chi);
    free(pRandWerte[i1].det_dChi);
  }
  free(pRandWerte);
}

void GenericAnsatzFunction::completeElementList(){
  unsigned int  k;                  // element under consideration

  // update element list
  for (unsigned int i=0; i<waveletList.sizeWaveletList; ++i) {
    for (unsigned int j=0; j<waveletList.W[i].noElements; ++j) {
      k = waveletList.W[i].element[j];
      if (elementTree.element[k].noWavelets%delta == 0) {
        elementTree.element[k].wavelet = (unsigned int*) realloc(elementTree.element[k].wavelet,(elementTree.element[k].noWavelets+delta)*sizeof(unsigned int));
      }
      elementTree.element[k].wavelet[elementTree.element[k].noWavelets++] = i;
    }
  }

  // optimize memory for element list
  for (unsigned int i=0; i<elementTree.totalSizeElementList; ++i)
    elementTree.element[i].wavelet = (unsigned int*) realloc(elementTree.element[i].wavelet,elementTree.element[i].noWavelets*sizeof(unsigned int));
  return;
}
///@todo implement tau, kappa, Phi times Phi, include memset
