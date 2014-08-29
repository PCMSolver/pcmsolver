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

#ifndef CPCMSOLVER_HPP
#define CPCMSOLVER_HPP

#include <iosfwd>
#include <string>
#include <vector>

#include "Config.hpp"

#include <Eigen/Dense>

class Cavity;

#include "IGreensFunction.hpp"
#include "PCMSolver.hpp"
#include "SolverData.hpp"
#include "SolverFactory.hpp"

/*! \file CPCMSolver.hpp
 *  \class CPCMSolver
 *  \brief Solver for conductor-like approximation: C-PCM (COSMO)
 *  \author Roberto Di Remigio
 *  \date 2013
 */

class CPCMSolver : public PCMSolver
{
private:
    bool builtIsotropicMatrix;
    bool builtAnisotropicMatrix;
    bool hermitivitize_;
    double correction_;
    Eigen::MatrixXd fullPCMMatrix;
    std::vector<Eigen::MatrixXd> blockPCMMatrix;
    void buildIsotropicMatrix(const Cavity & cavity);
    virtual std::ostream & printSolver(std::ostream & os);
public:
    CPCMSolver() {}
    CPCMSolver(IGreensFunction * gfInside, IGreensFunction * gfOutside, bool symm,
               double correction)
        : PCMSolver(gfInside, gfOutside), builtIsotropicMatrix(false),
          builtAnisotropicMatrix(false), hermitivitize_(symm), correction_(correction) {}
    virtual ~CPCMSolver() {}
    virtual void buildSystemMatrix(const Cavity & cavity);
    virtual void compCharge(const Eigen::VectorXd & potential, Eigen::VectorXd & charge,
                            int irrep = 0);
    void correction(double corr) {
        correction_ = corr;
    }
    friend std::ostream & operator<<(std::ostream & os, CPCMSolver & solver) {
        return solver.printSolver(os);
    }
};

namespace
{
    PCMSolver * createCPCMSolver(const solverData & _data)
    {
        return new CPCMSolver(_data.gfInside, _data.gfOutside, _data.hermitivitize,
                              _data.correction);
    }
    const std::string CPCMSOLVER("CPCM");
    const bool registeredCPCMSolver = SolverFactory::TheSolverFactory().registerSolver(
                                          CPCMSOLVER, createCPCMSolver);
}

#endif // CPCMSOLVER_HPP
