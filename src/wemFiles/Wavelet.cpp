#include "Wavelet.hpp"
#include "Constants.hpp"

void addElement(Wavelet* w, et_root E, double* newWeights,unsigned int noNewWeights, unsigned int index){
  unsigned int k, l;
  unsigned int todo = 0;
  for(k = 0; k < noNewWeights; ++k){
    if(newWeights[k] !=0 ) {
      todo = 1;
      break;
    }
  }
  if ( todo == 0 ) return;
  
  if (w->noElements != 0) for( k = 0; (k<w->noElements) && (w->element[k]!=index); ++k);

  if( k < w->noElements ){
    // element already in the support of this wavelet, add to weights
    for( l = 0; l < noNewWeights; ++l){
      w->weight[k*noNewWeights+l] += newWeights[l];
    }
  } else {

    if ( w->noElements%delta == 0 ){      
      w->element = (unsigned int*) realloc(w->element, (w->noElements+delta)*sizeof(unsigned int));
      w->weight = (double*) realloc (w->weight, (w->noElements+delta)*noNewWeights*sizeof(double));
#ifdef DEBUG
      // make sure that we do not have rubbish in our data
      for( l = 0; l < delta; ++l){
        w->element[w->noElements+l] = 0;
      }
      for( l = 0; l < delta*noNewWeights; ++l){
        w->weight[w->noElements*noNewWeights+l] = 0;
      }
#endif
    }
    for( l = 0; l < noNewWeights; ++l){
      w->weight[w->noElements*noNewWeights+l] = newWeights[l];
    }
    w->element[w->noElements++] = index;
  }
}
