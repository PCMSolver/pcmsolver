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

#ifndef PWLSOLVER_HPP
#define PWLSOLVER_HPP

#include <iosfwd>
#include <string>

#include "Config.hpp"

#include <Eigen/Dense>

extern "C"
{
//#include "vector3.h"
//#include "sparse2.h"
//#include "intvector_pwl.h"
//#include "basis_pwl.h"
}

#include "IGreensFunction.hpp"
#include "SolverData.hpp"
#include "SolverFactory.hpp"
#include "WEMSolver.hpp"

/*! \file PWLSolver.hpp
 *  \class PWLSolver
 *  \brief Wavelet solver, piecewise linear.
 *  \author Luca Frediani
 *  \date 2012
 */

class PWLSolver : public WEMSolver
{
private:
    virtual void initPointers();
    virtual std::ostream & printSolver(std::ostream & os);
public:
    PWLSolver(IGreensFunction * gfInside_,
              IGreensFunction * gfOutside_) : WEMSolver(gfInside_, gfOutside_, FirstKind) {
        initPointers();
    }
    PWLSolver(IGreensFunction * gfInside_, IGreensFunction * gfOutside_,
              int integralEquation_) : WEMSolver(gfInside_, gfOutside_, integralEquation_) {
        initPointers();
    }
    //PWLSolver(Section solver);
    virtual ~PWLSolver();
    friend std::ostream & operator<<(std::ostream & os, PWLSolver & solver) {
        return solver.printSolver(os);
    }
private:
    virtual void initInterpolation();
    virtual void constructWavelets();
    virtual void constructSi();
    virtual void constructSe();
    virtual void solveFirstKind(const Eigen::VectorXd & potential,
                                Eigen::VectorXd & charge);
    virtual void solveSecondKind(const Eigen::VectorXd & potential,
                                 Eigen::VectorXd & charge);
    virtual void solveFull(const Eigen::VectorXd & potential, Eigen::VectorXd & charge);
    element_pwl *elementTree; //*E_; Hierarchical element list
    wavelet_pwl *waveletList; //*W_; List of wavelets
};

namespace
{
    PCMSolver * createPWLSolver(const solverData & _data)
    {
        return new PWLSolver(_data.gfInside, _data.gfOutside, _data.integralEquation);
    }
    const std::string PWLSOLVER("LINEAR");
    const bool registeredPWLSolver = SolverFactory::TheSolverFactory().registerSolver(
                                         PWLSOLVER, createPWLSolver);
}

#endif // PWLSOLVER_HPP
