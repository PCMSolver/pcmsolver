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
 *     PCMSolver API, see: <https://repo.ctcc.no/projects/pcmsolver>
 */
/* pcmsolver_copyright_end */

#ifndef PCMSOLVER_HPP
#define PCMSOLVER_HPP

#include <iosfwd>

#include "Config.hpp"


class Cavity;

#include "IGreensFunction.hpp"

/*!
 * \file PCMSolver.hpp
 * \class PCMSolver
 * \brief Abstract Base Class for solvers inheritance hierarchy.
 * \author Luca Frediani
 * \date 2011
 */

class PCMSolver
{
protected:
    IGreensFunction * greenInside_;
    IGreensFunction * greenOutside_;
    bool allocated;
    virtual std::ostream & printSolver(std::ostream & os) = 0;
public:
    PCMSolver() {}
    PCMSolver(IGreensFunction * gfInside,
              IGreensFunction * gfOutside) : greenInside_(gfInside), greenOutside_(gfOutside),
        allocated(true) {}
    virtual ~PCMSolver() {
        if (allocated) {
            delete greenInside_;
            delete greenOutside_;
        }
    }

    IGreensFunction * greenInside() const {
        return greenInside_;
    }
    IGreensFunction * greenOutside() const {
        return greenOutside_;
    }

    /*! \brief Calculation of the PCM matrix.
     *  \param[in] cavity the cavity to be used.
     */
    virtual void buildSystemMatrix(const Cavity & cavity) = 0;
    /*! \brief Computation of ASC given the MEP.
    *  \param[in] potential the vector containing the MEP at cavity points.
     *  \param[out] charge the vector containing the ASC at cavity points.
     *  \param[in] irrep the irreducible representation of the MEP and ASC.
     *
     *  Given the MEP for a certain irrep, computes the corresponding ASC.
     *  By default, we expect the totally symmetric irrep to be needed,
     *  as in energy calculations..
     */
    virtual void compCharge(const Eigen::VectorXd & potential, Eigen::VectorXd & charge,
                            int irrep = 0) = 0;
    friend std::ostream & operator<<(std::ostream & os, PCMSolver & solver) {
        return solver.printSolver(os);
    }
};

#endif // PCMSOLVER_HPP
