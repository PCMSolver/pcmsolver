#ifndef WAVELETS_HPP
#define WAVELETS_HPP

/**
 * @file Wavelet.hpp
 *
 * @brief class declaration of wavelet
 */

#include "ElementTree.hpp"
#include <stdlib.h>

class Wavelet {
  public:
    unsigned int level;      ///< level of wavelet 
    unsigned int noElements; ///< elements in the support
    unsigned int noSons;     ///< number of sons, in the ConAF case we have always 4, in the LinAF case we might sometimes have only 2
    unsigned int *element;   ///< vector of indeces to elements
    unsigned int *son;       ///< vector of indeces to sons
    double *weight;          ///< the weight of the elements

    Wavelet():level(0), noElements(0), noSons(0),element(NULL),son(NULL),weight(NULL){};

};

class WaveletRoot{
  public:
    unsigned int sizeWaveletList; ///< the size of the waveletList
    Wavelet *W;                   ///< pointer to the wavelets

    WaveletRoot():sizeWaveletList(0),W(NULL){};
};

// adds one element with newWeights to the support of w
void addElement(Wavelet* w, et_root E, double* newWeights,unsigned int noWeights, unsigned int index);
#endif
