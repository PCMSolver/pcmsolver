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
#include <memory>
#include <string>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

class Cavity;

#include "Factory.hpp"
#include "IGreensFunction.hpp"
#include "PCMSolver.hpp"
#include "SolverData.hpp"

/*! \file CPCMSolver.hpp
 *  \class CPCMSolver
 *  \brief Solver for conductor-like approximation: C-PCM (COSMO)
 *  \author Roberto Di Remigio
 *  \date 2013
 */

class CPCMSolver : public PCMSolver
{
public:
    CPCMSolver() {}
    /*! \brief Construct solver from two shared_ptr to Green's functions
     *  \param[in] gf_i Green's function inside the cavity
     *  \param[in] gf_o Green's function outside the cavity
     *  \param[in] symm whether the system matrix has to be symmetrized
     *  \param[in] corr factor to correct the conductor results
     */
    CPCMSolver(const SharedIGreensFunction & gf_i, const SharedIGreensFunction & gf_o,
            bool symm, double corr)
        : PCMSolver(gf_i, gf_o), hermitivitize_(symm), correction_(corr) {}
    /*! \brief Construct solver from two raw pointers to Green's functions
     *  \param[in] gf_i Green's function inside the cavity
     *  \param[in] gf_o Green's function outside the cavity
     *  \param[in] symm whether the system matrix has to be symmetrized
     *  \param[in] corr factor to correct the conductor results
     *  \warning gf_i and gf_o will be deallocated automatically when the solver object goes out of scope,
     *  since they are wrapped in a std::shared_ptr
     */
    CPCMSolver(IGreensFunction * gf_i, IGreensFunction * gf_o, bool symm, double corr)
        : PCMSolver(gf_i, gf_o), hermitivitize_(symm), correction_(corr) {}
    virtual ~CPCMSolver() {}
    friend std::ostream & operator<<(std::ostream & os, CPCMSolver & solver) {
        return solver.printSolver(os);
    }
private:
    /*! Whether the system matrix has to be symmetrized */
    bool hermitivitize_;
    /*! Correction for the conductor results */
    double correction_;
    /*! PCM matrix, not symmetry blocked */
    Eigen::MatrixXd fullPCMMatrix_;
    /*! PCM matrix, symmetry blocked form */
    std::vector<Eigen::MatrixXd> blockPCMMatrix_;

    /*! \brief Calculation of the PCM matrix
     *  \param[in] cavity the cavity to be used
     */
    virtual void buildSystemMatrix_impl(const Cavity & cavity) override;
    /*! \brief Returns the ASC given the MEP and the desired irreducible representation
     *  \param[in] potential the vector containing the MEP at cavity points
     *  \param[in] irrep the irreducible representation of the MEP and ASC
     */
    virtual Eigen::VectorXd computeCharge_impl(const Eigen::VectorXd & potential,
            int irrep = 0) const override;
    virtual std::ostream & printSolver(std::ostream & os);
};

namespace
{
    PCMSolver * createCPCMSolver(const solverData & data, const SharedIGreensFunction & gf_i,
            const SharedIGreensFunction & gf_o)
    {
        return new CPCMSolver(gf_i, gf_o, data.hermitivitize, data.correction);
    }
    const std::string CPCMSOLVER("CPCM");
    const bool registeredCPCMSolver =
        Factory<PCMSolver, solverData, SharedIGreensFunction, SharedIGreensFunction>::TheFactory().registerObject(
                                         CPCMSOLVER, createCPCMSolver);
}

#endif // CPCMSOLVER_HPP
