#ifndef GENERICANSATZFUNCTION_HPP
#define GENERICANSATZFUNCTION_HPP

/**
 * @file GenericAnsatzFunction.hpp
 *
 * @brief class for a generic ansatz functions
 */
#include "BoundingBox.hpp"
#include "Constants.hpp"
#include "Cubature.hpp"
#include "ElementTree.hpp"
#include "Interpolation.hpp"
#include "Randwerte.hpp"
#include "SparseMatrix.hpp"
#include "Wavelet.hpp"

class Vector3;
class Vector2;

#define NEWTON 1
#define LAPLACE 2

class GenericAnsatzFunction {
public:

  Vector3 *nodeList;               ///< points for the single scale basis
  unsigned int **pElementList;       ///< elements in the single scale basis

  unsigned int nLevels;        ///< 2^nLevels*2^nLevels elements per patch
  unsigned int nFunctions;      ///< size of F, the element list of finest level
  unsigned int nNodes;        ///< size of P, the point list
  unsigned int nPatches;              ///< number of patches

  unsigned int noPhi;                ///< number of basis functions per element
  unsigned int quadratureLevel_;     ///< the degree for the RHS quadrature
  unsigned int td;                   ///< tilde_d parameter - the dual order parameter
  double dp;                         ///< d' parameter from paper, dp in (d, tilde_d+r) for postproc, or in (d, tilde_d+2*q) for compression

  Randwerte* pRandWerte;             ///< precalulated boundary values on gamma

  Interpolation *interCoeff;         ///< coefficients of surface interpolation

  et_root elementTree;               ///< root of hierarchical element tree

  WaveletRoot waveletList;           ///< list of wavelets

  BoundingBox *B;                    ///< the surrounding box

  SparseMatrix *G;                   ///< gram matrix - in the linear case

  /// constructor of the generic class AnsatzFunction
  GenericAnsatzFunction() {};

  /// generate point and patch list in C Format for the single scale basis and the ElementTree
  unsigned int genNet(Vector3 ***U);

  void completeElementList();
  virtual void simplifyWaveletList() = 0;
  
  int print_geometry(double * rho, char* dname);
  int printElement(char *dname, int element);

  /**
   * @brief reserves space for border-values that might be used again, however
   * the initialization is done in the function integrateNoProblem
   */
  void initRandwerte(int g_max);
  /**
   * @brief resets the space for border-values that might be used again
   */
  void resetRandwerte(int g_max);
  /**
   * @brief frees the space previously allocated for border values
   */
  void freeRandwerte();

  /**
   * @brief checks two elements with index e1 and e2 for a common edge (return
   * 3) or vertex (return 4), ind1 and ind2 are the corresponding rotations of
   * the elements. In case of no common intersection returns 1. In case of same
   * patch, the comparison is done via indeces. In case of different patches the
   * comparison is done via differences
   *
   * @param e1 e2 index of elements to compare
   * @param ind1 ind2 indeces for the rotations
   */
  unsigned int compare(unsigned int e1, unsigned int e2, unsigned int *ind1,
      unsigned int *ind2);
  double distance(int e1, int e2);


  void elementElementInteraction(double *c, unsigned int ind1, unsigned int ind2, 
      double alpha, Cubature *Q, double SingleLayer(Vector3, Vector3), 
      double DoubleLayer(Vector3, Vector3, Vector3), double Identity);

  virtual void quadratureGrade(signed int *g1, signed int*g2, int level1, 
      int level2, double dist, double alpha) = 0;  

  //integrate functions
  /// integration function for - no problem quadrature - modified scalar product
  virtual void integrateNoProblem(double *c, unsigned int i1, unsigned int i2, unsigned int g1, Cubature *Q1, Cubature *Q2, double SingleLayer(Vector3, Vector3), double DoubleLayer(Vector3, Vector3, Vector3)) = 0;

  /// integration function for - same patches
  virtual void integratePatch(double *c, unsigned int i1, Cubature * Q, double SingleLayer(Vector3, Vector3), double DoubleLayer(Vector3, Vector3, Vector3), double Identity) = 0;

  /// integration function for - common edge
  virtual void integrateEdge(double *c, unsigned int i1, unsigned int i2, unsigned int ind_s, unsigned int ind_t, Cubature *Q, double SingleLayer(Vector3, Vector3), double DoubleLayer(Vector3, Vector3, Vector3), double Identity) = 0;

  /// integration function for - common node in origin
  virtual void integratePoint(double *c, unsigned int i1, unsigned int i2, unsigned int ind_s,
		unsigned int ind_t, Cubature *Q, double SingleLayer(Vector3, Vector3), double DoubleLayer(Vector3, Vector3, Vector3), double Identity) = 0;

  signed int searchIntegral(intvector *I, unsigned int i);
  void setIntegral(intvector *I, unsigned int i, double* z);
  virtual void calculateIntegral(double *a, double *c) = 0;

  virtual void permutate (double *a, double *b) = 0;

  virtual void calculateYRHS(double **y, int i) = 0;
  virtual void calculateCRHS(double *c, double w, Vector2 xi) = 0;
  virtual double calculateUEnergy(double *u, Vector2 xi, unsigned int zi) = 0;

  /// construction of WaveletList
  virtual void generateWaveletList() = 0;

  virtual void setQuadratureLevel() = 0;

  unsigned int compression(SparseMatrix *T);
  unsigned int postProc(SparseMatrix *T);

  virtual void computeBoundingBoxes() = 0;
  virtual void freeBoundingBoxes() = 0;

  virtual unsigned int waveletWaveletCriterion(unsigned int ind1, unsigned int ind2, double c1, double c2) = 0;

  virtual void createGram(unsigned int size, unsigned int maxRowNum) = 0;
  virtual void freeGram() = 0;

  void inv_A_times_x(double* v);

  void precond(double* a, double* b);

  virtual void tdwt(double *v) = 0;
  virtual void  dwt(double *v) = 0;
  
  int freeElementTree();
  ~GenericAnsatzFunction(){};
};
#endif
