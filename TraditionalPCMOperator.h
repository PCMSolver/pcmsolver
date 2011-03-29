#ifndef TRADITIONALPCMOPERATOR_H
#define TRADITIONALPCMOPERATOR_H


//#include "BEMInterface.h"
#include <string>
#include <vector>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/LU>
//#include "GreensFunctionInterface.h"


// #include "GreensFunctionInterface.h"


USING_PART_OF_NAMESPACE_EIGEN;

//class TraditionalPCMOperator:public BEMInterface{
class TraditionalPCMOperator{


public:

  TraditionalPCMOperator(double permitivity){
    PCM_E = permitivity;
  }

  virtual ~TraditionalPCMOperator() {};

  /* @brief Clone the interface type
   *
   * @return newly created BEMInterface-object
   *
   */
  //  virtual BEMInterface *newBEMInterface();

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
  
  virtual void getCavityDefs(std::vector<Vector3d> &points_,
			     std::vector<Vector3d> &normals_,
			     std::vector<double> &areas_);
  

  /* @brief Constructs the system matrix according to some interface.
   *
   * Remember to call this method before computing charges etc. NOT
   * thread safe, i.e., only one system matrix per object.
   *
   */
  //  virtual void constructSystemMatrix(GreensFunctionInterface &iface);
  virtual void constructSystemMatrix();  

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
  //   virtual void solveForPotential();
  

  /* @brief Output
   *
   */
  virtual void printInfo(std::ostream &out);

  //  void test();

protected:

  double PCM_E;  
  std::vector<Vector3d> points_;
  std::vector<double> areas_;
  std::vector<Vector3d> correspondingSpheres_;
  std::vector<Vector3d> atoms_;
  std::vector<double> sphereRadii_;
  std::vector<Vector3d> normals_;
  MatrixXd systemMatrix_;
};

#endif
