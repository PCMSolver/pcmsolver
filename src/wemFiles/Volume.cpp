#include "GenericAnsatzFunction.hpp"
#include "Cubature.hpp"
#include "GaussSquare.hpp"

#include "stdio.h"
#include "string.h"

double calculateVolume(GenericAnsatzFunction *af){
  unsigned int n = 1<<af->nLevels;
  double h = 1./n;
  double V = 0.0;

  Cubature *Q;

  Vector3 tPoint, ntPoint;
  Vector2 t; 
  initGaussSquare(&Q, af->quadratureLevel_+1);

  for(unsigned int i = af->nPatches*(n*n-1)/3; i < af->elementTree.totalSizeElementList; ++i){
    for(unsigned int k = 0; k < Q[af->quadratureLevel_].noP; ++k){
      t.x = h*(af->elementTree.element[i].index_s+Q[af->quadratureLevel_].xi[k].x);
      t.y = h*(af->elementTree.element[i].index_t+Q[af->quadratureLevel_].xi[k].y);
      tPoint = af->interCoeff->Chi(t,af->elementTree.element[i].patch);
      ntPoint = af->interCoeff->n_Chi(t,af->elementTree.element[i].patch);
      
      V += Q[af->quadratureLevel_].weight[k]*vector3Dot(tPoint, ntPoint);
    }
  }
  freeGaussSquare(&Q, af->quadratureLevel_+1);
  return (h*h*V/3);
}
