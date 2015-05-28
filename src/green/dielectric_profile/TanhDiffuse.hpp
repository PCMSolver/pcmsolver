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

#ifndef TANHDIFFUSE_HPP
#define TANHDIFFUSE_HPP

#include <iosfwd>
#include <tuple>

#include "Config.hpp"

/*! \file TanhDiffuse.hpp
 *  \class TanhDiffuse
 *  \brief A tanh dielectric profile as in \cite Frediani2004a
 *  \author Roberto Di Remigio
 *  \date 2014
 *  \note width_ is the parameter from user input. The division by 6.0 is
 *  to keep consistency with \cite Frediani2004a
 */

class TanhDiffuse
{
private:
    /// Dielectric constant on the left of the interface
    double epsilon1_;
    /// Dielectric constant one the right of the interface
    double epsilon2_;
    /// Width of the transition layer
    double width_;
    /// Center of the transition layer
    double center_;
    /*! Returns value of dielectric profile at given point
     *  \param[in] point where to evaluate the profile
     */
    double value(double point) const {
	double effWidth = width_ / 6.0;
        double epsPlus = (epsilon1_ + epsilon2_) / 2.0;
        double epsMinus = (epsilon2_ - epsilon1_) / 2.0;
        double tanh_r = std::tanh((point - center_) / effWidth);
        return (epsPlus + epsMinus * tanh_r); //epsilon(r)
    }
    /*! Returns value of derivative of dielectric profile at given point
     *  \param[in] point where to evaluate the derivative
     */
    double derivative(double point) const {
	double effWidth = width_ / 6.0;
        double factor = (epsilon1_ - epsilon2_) / (2.0 * effWidth);
        double tanh_r = std::tanh((point - center_) / effWidth);
        return (factor * (1 - std::pow(tanh_r, 2))); //first derivative of epsilon(r)
    }
    std::ostream & printObject(std::ostream & os)
    {
        os << "Permittivity inside  = " << epsilon1_ << std::endl;
        os << "Permittivity outside = " << epsilon2_ << std::endl;
        os << "Profile width        = " << width_ / 6.0    << " AU" << std::endl;
        os << "Profile center       = " << center_   << " AU";
        return os;
    }
public:
    TanhDiffuse() {}
    TanhDiffuse(double e1, double e2, double w, double c) :
        epsilon1_(e1), epsilon2_(e2), width_(w), center_(c) {}
    /*! Returns a tuple holding the permittivity and its derivative
     *  \param[in]   r evaluation point
     */
    std::tuple<double, double> operator()(const double r) const
    {
        return std::make_tuple(value(r), derivative(r));
    }
    double epsilon1() const { return epsilon1_; }
    double epsilon2() const { return epsilon2_; }
    double width() const { return (width_ / 6.0); }
    double center() const { return center_; }
    friend std::ostream & operator<<(std::ostream & os, TanhDiffuse & th) {
        return th.printObject(os);
    }
};

#endif // TANHDIFFUSE_HPP