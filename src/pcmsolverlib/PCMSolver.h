/*! \file PCMSolver.h 
\brief PCM solver
*/


#ifndef PCMSOLVER
#define PCMSOLVER

#include <string>
#include <vector>
#include <iostream>
#include <complex>

using namespace std;

#include "GreensFunction.h"
#include "UniformDielectric.h"
#include "GreensFunctionSum.h"

#include "MetalSphere.h"

template <class Green>
class PCMSolver{
 public:
    PCMSolver(Green &greensFunction_){greensFunction = &greensFunction_;};
    ~PCMSolver(){};
    Green &getGreensFunction();
 private:
    Green *greensFunction;
};
/*
class PCMSolver{
 public:
    PCMSolver(GreensFunction &GF_){GF = &GF_;};
    ~PCMSolver(){};
    GreensFunction &getGreensFunction();
 private:
    GreensFunction *GF;
*/
    /* 
protected:
    double PCM_E;  
    std::vector<Vector3d> points_;
    std::vector<double> areas_;
    std::vector<Vector3d> correspondingSpheres_;
    std::vector<Vector3d> atoms_;
    std::vector<double> sphereRadii_;
    std::vector<Vector3d> normals_;
    MatrixXd systemMatrix_;
    */
#endif
