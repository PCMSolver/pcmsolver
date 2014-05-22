#ifndef CONANSATZFUNCTION_HPP
#define CONANSATZFUNCTION_HPP

#include "GenericAnsatzFunction.hpp"

/**
 * @file ConAnsatzFunction.hpp
 *
 * @brief class for a constant Ansatz functions
 * @note contains all the functions that are specific to a Constant Ansatz Function
 */
class ConAnsatzFunction : public GenericAnsatzFunction {
  public:

  // characteristic variable, found only in this function
  BoundingBoxSquare *B2;
  
  /// constructor of the class ConAnsatzFunction
  ConAnsatzFunction();
  ConAnsatzFunction(unsigned int p, unsigned int m, unsigned int nf, Vector3*** pPointsIn);

  void quadratureGrade(signed int *g1, signed int*g2, int level1, int level2, double dist, double alpha);

  //integrate functions

  // integration function for - no problem quadrature - modified scalar product
  void integrateNoProblem(double *c, unsigned int i1, unsigned int i2, unsigned int g1, Cubature *Q1, Cubature *Q, double SingleLayer(Vector3, Vector3), double DoubleLayer(Vector3, Vector3, Vector3));

  // integration function for - same patches
  void integratePatch(double *c, unsigned int i1, Cubature * Q, double SingleLayer(Vector3, Vector3), double DoubleLayer(Vector3, Vector3, Vector3), double Identity);

  // integration function for - common edge
  void integrateEdge(double *c, unsigned int i1, unsigned int i2, unsigned int ind_s, unsigned int ind_t, Cubature *Q, double SingleLayer(Vector3, Vector3), double DoubleLayer(Vector3, Vector3, Vector3), double Identity);

  // integration function for - common node in origin
  void integratePoint(double *c, unsigned int i1, unsigned int i2, unsigned int ind_s,
		unsigned int ind_t, Cubature *Q, double SingleLayer(Vector3, Vector3), double DoubleLayer(Vector3, Vector3, Vector3), double Identity);
  
  void calculateIntegral(double *a, double *c);
  
  void calculateYRHS(double **y, int i);
  void calculateCRHS(double *c, double w, Vector2 xi);
  double calculateUEnergy(double *u, Vector2 xi, unsigned int zi);

  void permutate(double* a, double* b);
  
  // construction of WaveletList
  void generateWaveletList();
  
  void setQuadratureLevel();

  void simplifyWaveletList(); 

  void computeBoundingBoxes();
  void freeBoundingBoxes();
  
  unsigned int waveletWaveletCriterion(unsigned int ind1, unsigned int ind2, double c1, double c2);

  void createGram(unsigned int size, unsigned int maxRowNum);
  void freeGram();
  
  void tdwt(double *v);
  void  dwt(double *v);

  ~ConAnsatzFunction();
};
#endif
