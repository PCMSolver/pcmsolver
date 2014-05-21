#ifndef _ELEMENT_TREE_HPP
#define _ELEMENT_TREE_HPP

class Vector3;

typedef struct {
  double **value;              ///< the value of the integral
  unsigned int *index;         ///< the index
  unsigned int integralNumber; ///< the number of integrals stored
} intvector;


typedef struct _et_node {
  Vector3        midpoint;   ///< midpoint of circumscribed circle
  int            father;     ///< father element
  int            son[4];     ///< son elements
  _et_node       *up;         ///< pointer to the upper neighbour
  _et_node       *down;       ///< pointer to the bottom neighbour
  _et_node       *left;       ///< pointer to the left neighbour
  _et_node       *right;      ///< pointer to the right neighbour
  double         radius;     ///< radius of the circumscribed circle
  unsigned int   level;      ///< level of element
  int            number;     ///< index of element on current level
  unsigned int   patch;      ///< index of the patch of the element
  unsigned int   index_s;    ///< index of element in s direction
  unsigned int   index_t;    ///< index of element in t direction
  unsigned int   vertex[4];  ///< indeces of verteces
  unsigned int  *wavelet;    ///< indeces of the wavelets that have this element in their support
  unsigned int   noWavelets; ///< number of wavelets that have this element in their support
  intvector	interaction;     ///< sparse vector representation for integrals
} et_node;

typedef struct {
  unsigned int totalSizeElementList; ///< the total size of the element list
  int     nop;                       ///< the number of patches
  et_node *element;                  ///< pointer to the element vector
} et_root;

#endif

