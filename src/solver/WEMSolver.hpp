/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013 Roberto Di Remigio, Luca Frediani and contributors
 *     
 *     This file is part of PCMSolver.
 *     
 *     PCMSolver is free software: you can redistribute it and/or modify       
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *     
 *     PCMSolver is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *     
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *     
 *     For information on the complete list of contributors to the
 *     PCMSolver API, see: <http://pcmsolver.github.io/pcmsolver-doc>
 */
/* pcmsolver_copyright_end */

#ifndef WEMSOLVER_HPP
#define WEMSOLVER_HPP

#include <iosfwd>
#include <string>

#include "Config.hpp"

#include <Eigen/Dense>

#include "SparseMatrix.hpp"
#include "GenericAnsatzFunction.hpp"
#include "LinAnsatzFunction.hpp"
#include "ConAnsatzFunction.hpp"

class Vector2;
class Vector3;
class Cavity;
class WaveletCavity;

#include "IGreensFunction.hpp"
#include "SolverData.hpp"
#include "SolverFactory.hpp"
#include "PCMSolver.hpp"

/*! \file WEMSolver.hpp
 *  \class WEMSolver
 *  \brief WEMSolver
 *  \author Luca Frediani and Monica Bugeanu
 *  \date 2014
 */

class WEMSolver : public PCMSolver
{
private:
    void initWEMMembers();
    virtual std::ostream & printSolver(std::ostream & os);
public:
    WEMSolver(IGreensFunction * gfInside_, IGreensFunction * gfOutside_,
              std::string modelType, int integralEquation_ = SecondKind )
        : PCMSolver(gfInside_, gfOutside_), integralEquation(integralEquation_) {
        if(modelType =="Linear"){
          af = new LinAnsatzFunction();
          interpolationGrade = 2;
          interpolationType = 1;
        }
        else if(modelType == "Wavelet"){
          af = new ConAnsatzFunction();
          interpolationGrade = 1;
          interpolationType = 1;
        }
        initWEMMembers();
    }
//                WEMSolver(const Section & solver);
    virtual ~WEMSolver();
     Interpolation* getT_() {
        return af->interCoeff;
    }
    int getQuadratureLevel() {
        return af->quadratureLevel_;
    }
    virtual void buildSystemMatrix(const Cavity & cavity);
    virtual void compCharge(const Eigen::VectorXd & potential, Eigen::VectorXd & charge,
                            int irrep = 0);
    friend std::ostream & operator<<(std::ostream & os, WEMSolver & solver) {
        return solver.printSolver(os);
    }
protected:
    void constructSystemMatrix();
    void uploadCavity(const WaveletCavity &
                              cavity); // different interpolation
    virtual void constructSi();
    virtual void constructSe();
    virtual void solveFirstKind(const Eigen::VectorXd & potential,
                                Eigen::VectorXd & charge);
    virtual void solveSecondKind(const Eigen::VectorXd & potential,
                                 Eigen::VectorXd & charge);
    virtual void solveFull(const Eigen::VectorXd & potential,
                           Eigen::VectorXd & charge);
    virtual void constructWavelets();
    virtual void initInterpolation();

    double threshold;
    SparseMatrix S_i_, S_e_; // System matrices
    bool systemMatricesInitialized_;
    GenericAnsatzFunction *af;

//    unsigned int quadratureLevel_;
    Vector3 *** pointList; // the old U
//    Vector3 *nodeList; //*P_; --     // Point list
//    unsigned int **elementList; //**F_;     // Element list
//    unsigned int nNodes; //np_; --    // Number of knot points or something
//    unsigned int nFunctions; //nf_; --    // Number of ansatz functions
//    unsigned int nPatches; // p_; --    // Number of points
//    unsigned int nLevels; //M_; --    // Patch level (2**M * 2**M elements per patch)
//    int nQuadPoints; // nPoints_;    // Number of quadrature points
    
    double apriori1_, aposteriori1_;    // System matrix sparsities
    double apriori2_, aposteriori2_;    // System matrix sparsities
    int integralEquation;

    enum EquationType {FirstKind, SecondKind, Full};
};

namespace
{
    PCMSolver * createWEMSolver(const solverData & _data)
    {
        return new WEMSolver(_data.gfInside, _data.gfOutside, _data.modelType, _data.integralEquation);
    }
    const std::string WEMSOLVER("Wavelet");
    const bool registeredWEMSolver = SolverFactory::TheSolverFactory().registerSolver(
        WEMSOLVER,createWEMSolver);
}

#endif // WEMSOLVER_HPP
