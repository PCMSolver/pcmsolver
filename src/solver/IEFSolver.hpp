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

#include <Eigen/Core>

class Cavity;
class IGreensFunction;

#include "Factory.hpp"
#include "PCMSolver.hpp"
#include "SolverData.hpp"

/*! \file IEFSolver.hpp
 *  \class IEFSolver
 *  \brief IEFPCM, collocation-based solver
 *  \author Luca Frediani, Roberto Di Remigio
 *  \date 2011, 2015
 */

class IEFSolver : public PCMSolver
{
public:
    IEFSolver() {}
    /*! \brief Construct solver from two shared_ptr to Green's functions
     *  \param[in] symm whether the system matrix has to be symmetrized
     */
    IEFSolver(bool symm) : PCMSolver(), hermitivitize_(symm) {}
    virtual ~IEFSolver() {}
    /*! \brief Builds PCM matrix for an anisotropic environment
     *  \param[in] cavity the cavity to be used.
     *  \param[in] gf_i Green's function inside the cavity
     *  \param[in] gf_o Green's function outside the cavity
     */
    void buildAnisotropicMatrix(const Cavity & cavity, const IGreensFunction & gf_i, const IGreensFunction & gf_o);
    /*! \brief Builds PCM matrix for an isotropic environment
     *  \param[in] cavity the cavity to be used.
     *  \param[in] gf_i Green's function inside the cavity
     *  \param[in] gf_o Green's function outside the cavity
     */
    void buildIsotropicMatrix(const Cavity & cavity, const IGreensFunction & gf_i, const IGreensFunction & gf_o);
    friend std::ostream & operator<<(std::ostream & os, IEFSolver & solver) {
        return solver.printSolver(os);
    }
private:
    /*! Whether the system matrix has to be symmetrized */
    bool hermitivitize_;
    /*! PCM matrix, not symmetry blocked */
    Eigen::MatrixXd fullPCMMatrix_;
    /*! PCM matrix, symmetry blocked form */
    std::vector<Eigen::MatrixXd> blockPCMMatrix_;

    /*! \brief Calculation of the PCM matrix
     *  \param[in] cavity the cavity to be used
     *  \param[in] gf_i Green's function inside the cavity
     *  \param[in] gf_o Green's function outside the cavity
     */
    virtual void buildSystemMatrix_impl(const Cavity & cavity, const IGreensFunction & gf_i, const IGreensFunction & gf_o) __override;
    /*! \brief Returns the ASC given the MEP and the desired irreducible representation
     *  \param[in] potential the vector containing the MEP at cavity points
     *  \param[in] irrep the irreducible representation of the MEP and ASC
     */
    virtual Eigen::VectorXd computeCharge_impl(const Eigen::VectorXd & potential,
            int irrep = 0) const __override;
    virtual std::ostream & printSolver(std::ostream & os) __override;
};

namespace
{
    PCMSolver * createIEFSolver(const solverData & data)
    {
        return new IEFSolver(data.hermitivitize);
    }
    const std::string IEFSOLVER("IEFPCM");
    const bool registeredIEFSolver =
        Factory<PCMSolver, solverData>::TheFactory().registerObject(IEFSOLVER, createIEFSolver);
}

#endif // IEFSOLVER_HPP
