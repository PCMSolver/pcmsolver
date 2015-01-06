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

#ifndef IEFSOLVER_HPP
#define IEFSOLVER_HPP

#include <iosfwd>
#include <string>
#include <vector>

#include "Config.hpp"


class Cavity;

#include "IGreensFunction.hpp"
#include "PCMSolver.hpp"
#include "SolverData.hpp"
#include "SolverFactory.hpp"

/*!
 * \file IEFSolver.hpp
 * \class IEFSolver
 * \brief Traditional solver.
 * \author Luca Frediani
 * \date 2011
 */

class IEFSolver : public PCMSolver
{
private:
    bool builtIsotropicMatrix;
    bool builtAnisotropicMatrix;
    bool hermitivitize_;
    Eigen::MatrixXd fullPCMMatrix;
    std::vector<Eigen::MatrixXd> blockPCMMatrix;
    /*! \brief Builds PCM matrix for an anisotropic environment
     *  \param[in] cavity the cavity to be used.
     */
    void buildAnisotropicMatrix(const Cavity & cavity);
    /*! \brief Builds PCM matrix for an isotropic environment
     *  \param[in] cavity the cavity to be used.
     */
    void buildIsotropicMatrix(const Cavity & cavity);
    virtual std::ostream & printSolver(std::ostream & os);
public:
    IEFSolver() {}
    IEFSolver(IGreensFunction * gfInside, IGreensFunction * gfOutside, bool symm)
        : PCMSolver(gfInside, gfOutside), builtIsotropicMatrix(false),
          builtAnisotropicMatrix(false), hermitivitize_(symm) {}
    virtual ~IEFSolver() {}
    virtual void buildSystemMatrix(const Cavity & cavity);
    virtual void computeCharge(const Eigen::VectorXd &potential, Eigen::VectorXd &charge,
            int irrep = 0);
    friend std::ostream & operator<<(std::ostream & os, IEFSolver & solver) {
        return solver.printSolver(os);
    }
};

namespace
{
    PCMSolver * createIEFSolver(const solverData & _data)
    {
        return new IEFSolver(_data.gfInside, _data.gfOutside, _data.hermitivitize);
    }
    const std::string IEFSOLVER("IEFPCM");
    const bool registeredIEFSolver = SolverFactory::TheSolverFactory().registerSolver(
                                         IEFSOLVER, createIEFSolver);
}

#endif // IEFSOLVER_HPP
