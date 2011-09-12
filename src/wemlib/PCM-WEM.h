/*
** PCM-WEM.h
** 
** Made by Ville Weijo
** Login   <weijo@weijo-desktop>
** 
** Started on  Wed Dec  9 16:27:19 2009 Ville Weijo
** Last update Wed Dec  9 16:27:19 2009 Ville Weijo
*/

#ifndef   	PCMWEM_H_
#define   	PCMWEM_H_

// All that is included...
extern "C"{
#include "vector3.h"
#include "sparse2.h"
#include "intvector.h"
#include "basis.h"
}
#include <vector>
#include <string>
#include <fstream>
#include <Eigen/Core>

// From the PCM module
#include "BEMInterface.h"
#include "GreensFunctionInterface.h"

using namespace Eigen;

// Class interface for the wavelet solver

class WaveletBEM: public BEMInterface{

 public:

  WaveletBEM();
  virtual ~WaveletBEM();

  // The following is from the BEMInterface

  /* @brief Clone the interface type
   *
   * @return newly created BEMInterface-object
   *
   */
  virtual BEMInterface *newBEMInterface();
  
  /* @brief Initializes the cavity
   *
   * Initializes cavity definitions from a file. This should be done
   * before any other operations.
   *
   * @param filename cavity file
   * @return false if everything went smoothly
   *
   */
  virtual bool readCavity(std::string &filename);


  /* @brief Returns essential cavity information. 
   *
   * @param points Quadrature points on the cavity
   * @param normals Normal of the cavity at a corresponding quadrature
   * point
   * @param areas Area/weight corresponding to a quadrature point
   *
   */
  virtual void getCavityDefs(std::vector<Vector3d> &points,
			     std::vector<Vector3d> &normals,
			     std::vector<double> &areas);

  /* @brief Constructs the system matrix according to some interface.
   *
   * Remember to call this method before computing charges etc. NOT
   * thread safe, i.e., only one system matrix per object.
   *
   */
  virtual void constructSystemMatrix(GreensFunctionInterface &iface);
  

  /* @brief Solve polarization charges according to some potential.
   *
   * Uses system matrix alread computed with constructSystemMatrix to
   * obtain polarization charges on the surface. Charges are returned
   * at quadrature points and they can be used to calculate
   * interaction energy directly by taking a dot product with the
   * potential.
   *
   * Thread safe.
   *
   * @param potential Potential calculated at the quadrature points
   * @param charges Resulting charges on surface, should have enough
   * space reserved to hold the numbers
   *
   */
  virtual void solveForPotential(VectorXd &potential, VectorXd &charges);
  

  /* @brief Writes charges to disk. 
   *
   * Writing is done in this class because one has to take into
   * account the order of quadrature etc. At the moment, charge on a
   * patch is averaged from the charges calculated at quadrature
   * points.
   *
   * @param Charges to be written to disk
   *
   */   
  void writeChargesToDisk(VectorXd &charges);


  /* @brief Output
   *
   */
  virtual void printInfo(std::ostream &out);

protected:

  // Parameters
  double eps_;
  unsigned int quadratureLevel_;



  // System matrices
  sparse2 S_i_, S_e_;
  bool systemMatricesInitialized_;
  // Point list
  vector3 *P_;
  // Element list
  unsigned int **F_;
  // Something 1
  vector3 ****T_;
  // Number of knot points or something
  unsigned int np_;
  // Number of ansatz functions
  unsigned int nf_;
  // Number of points 
  unsigned int p_;
  // Patch level (2**M * 2**M elements per patch)
  unsigned int M_;
  // Hierarchical element list
  element *E_;
  // List of wavelets
  wavelet *W_;
  // Number of quadrature points
  int nPoints_;

  // System matrix sparsities
  double apriori1_, aposteriori1_;
  double apriori2_, aposteriori2_;

};



#endif 	    /* !PCM-WEM_H_ */
