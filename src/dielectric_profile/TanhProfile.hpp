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

#ifndef TANHPROFILE_HPP
#define TANHPROFILE_HPP

#include <iosfwd>
#include <string>

#include "Config.hpp"

#include <Eigen/Dense>

/*! \file TanhProfile.hpp
 *  \class TanhProfile
 *  \brief A tanh dielectric profile as in \cite Frediani2004a
 *  \author Roberto Di Remigio
 *  \date 2014
 */

class TanhProfile
{
private:
    double epsilonLeft_;
    double epsilonRight_;
    double width_;
    double center_;
public:
    TanhProfile(double eL, double eR, double w, double c) :
        epsilonLeft_(eL), epsilonRight_(eR), width_(w), center_(c) {}
    double value(double point) const {
        double tanh_r = std::tanh((point - center_) / width_);
        return (0.5 * (epsilonLeft_ + epsilonRight_)
                + 0.5 * (epsilonRight_ - epsilonLeft_) * tanh_r);     //epsilon(r)
    }
    double derivative(double point) const {
        double tanh_r = std::tanh((point - center_) / width_);
        return (0.5 * (epsilonLeft_ - epsilonRight_)
                * ( 1 - tanh_r * tanh_r) / width_); //first derivative of epsilon(r)
    }
    /*! The permittivity profile of the transition layer
     *  \param[out]  e the value of the dielectric constant at point r
     *  \param[out] de the value of the derivative of the dielectric constant
     *                 at point r
     *  \param[in]   r evaluation point
     */
    void operator()(double & e, double & de, const double r) const {
        e  = value(r);
        de = derivative(r);
    }
};

#endif // TANHPROFILE_HPP
