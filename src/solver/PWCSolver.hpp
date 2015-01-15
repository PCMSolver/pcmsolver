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

#ifndef PWCSOLVER_HPP
#define PWCSOLVER_HPP

#include <iosfwd>
#include <string>

#include "Config.hpp"

#include <Eigen/Dense>

#include "SparseMatrix.hpp"
#include "GenericAnsatzFunction.hpp"
#include "ConAnsatzFunction.hpp"
#include "readPoints.hpp"

class Vector2;
class Vector3;
class Cavity;
class WaveletCavity;
struct Compression;

#include "IGreensFunction.hpp"
#include "SolverData.hpp"
#include "SolverFactory.hpp"
#include "PCMSolver.hpp"

/*! \file PWCSolver.hpp
 *  \class PWCSolver
 *  \brief Class describing a wavelet solver with piecewise constant wavelet functions 
 *  \author Luca Frediani and Monica Bugeanu
 *  \date 2014
 */

class PWCSolver : public PCMSolver
{
private:
    unsigned int interpolationGrade;
    unsigned int interpolationType;
    void initWEMMembers();
    virtual std::ostream & printSolver(std::ostream & os);
public:
    PWCSolver(IGreensFunction * gfInside_, IGreensFunction * gfOutside_, int integralEquation_ = SecondKind)
        : PCMSolver(gfInside_, gfOutside_), interpolationGrade(2), interpolationType(1), 
	af( new ConAnsatzFunction() ), integralEquation(integralEquation_) {
        initWEMMembers();
    }
    PWCSolver(IGreensFunction * gfInside_, IGreensFunction * gfOutside_, const Compression & _comp, int integralEquation_ = SecondKind)
        : PCMSolver(gfInside_, gfOutside_), interpolationGrade(2), interpolationType(1), 
	af( new ConAnsatzFunction(_comp) ), integralEquation(integralEquation_) {
        initWEMMembers();
    }
    virtual ~PWCSolver();
    Interpolation* getT_() {
        return af->interCoeff;
    }
    int getQuadratureLevel() {
        return af->quadratureLevel_;
    }
    virtual void buildSystemMatrix(const Cavity & cavity);
    virtual void computeCharge(const Eigen::VectorXd & potential, Eigen::VectorXd & charge,
                            int irrep = 0);
    friend std::ostream & operator<<(std::ostream & os, PWCSolver & solver) {
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
    GenericAnsatzFunction * af;

    Vector3 *** pointList; // the old U
    
    double apriori1_, aposteriori1_;    // System matrix sparsities
    double apriori2_, aposteriori2_;    // System matrix sparsities
    int integralEquation;

    enum EquationType {FirstKind, SecondKind, Full};
};

namespace
{
    PCMSolver * createPWCSolver(const solverData & _data)
    {
        return new PWCSolver(_data.gfInside, _data.gfOutside, _data.integralEquation);
    }
    const std::string PWCSOLVER("WAVELET"); // Stands for piecewise constant functions
    const bool registeredPWCSolver = SolverFactory::TheSolverFactory().registerSolver(
        PWCSOLVER, createPWCSolver);
}

#endif // PWCSOLVER_HPP
