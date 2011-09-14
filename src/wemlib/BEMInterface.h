#ifndef BEMINTERFACE_H
#define BEMINTERFACE_H

// Interface for a boundary element solver, e.g., wavelet solver

#include <string>
#include <vector>
#include <iostream>
#include <Eigen/Core>

#include "GreensFunctionInterface.h"

using namespace Eigen;

class BEMInterface{

 public:

  virtual ~BEMInterface() {};

  /* @brief Clone the interface type
   *
   * @return newly created BEMInterface-object
   *
   */
  virtual BEMInterface *newBEMInterface() = 0;

  /* @brief Initializes the cavity
   *
   * Initializes cavity definitions from a file. This should be done
   * before any other operations.
   *
   * @param filename cavity file
   * @return false if everything went smoothly
   *
   */
  virtual bool readCavity(std::string &filename) = 0;


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
			     std::vector<double> &areas) = 0;

  /* @brief Constructs the system matrix according to some interface.
   *
   * Remember to call this method before computing charges etc. NOT
   * thread safe, i.e., only one system matrix per object.
   *
   */
  virtual void constructSystemMatrix(GreensFunctionInterface &iface) = 0;
  

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
  virtual void solveForPotential(VectorXd &potential, VectorXd &charges) = 0;
  

  /* @brief Output
   *
   */
  virtual void printInfo(std::ostream &out) = 0;

  friend std::ostream& operator<<(std::ostream &out, BEMInterface &bem){
    bem.printInfo(out);
    return out;
  }

};

#endif
