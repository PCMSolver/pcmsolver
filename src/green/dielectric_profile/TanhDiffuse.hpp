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

#include "Config.hpp"

/*! \file TanhDiffuse.hpp
 *  \class TanhDiffuse
 *  \brief A tanh dielectric profile as in \cite Frediani2004a
 *  \author Roberto Di Remigio
 *  \date 2014
 */

class TanhDiffuse
{
private:
    double epsilon1_;
    double epsilon2_;
    double width_;
    double center_;
    double value(double point) const {
        double epsPlus = (epsilon1_ + epsilon2_) / 2.0;
        double epsMinus = (epsilon2_ - epsilon1_) / 2.0;
        double tanh_r = std::tanh((point - center_) / width_);
        return (epsPlus + epsMinus * tanh_r); //epsilon(r)
    }
    double derivative(double point) const {
        double factor = (epsilon1_ - epsilon2_) / (2.0 * width_);
        double tanh_r = std::tanh((point - center_) / width_);
        return (factor * (1 - std::pow(tanh_r, 2))); //first derivative of epsilon(r)
    }
    std::ostream & printObject(std::ostream & os)
    {
        os << "Permittivity inside  = " << epsilon1_ << std::endl;
        os << "Permittivity outside = " << epsilon2_ << std::endl;
        os << "Profile width        = " << width_    << std::endl;
        os << "Profile center       = " << center_;
        return os;
    }
public:
    TanhDiffuse() {}
    TanhDiffuse(double e1, double e2, double w, double c) :
        epsilon1_(e1), epsilon2_(e2), width_(w), center_(c) {}
    /*! The permittivity profile of the transition layer
     *  \param[out]  e the value of the dielectric constant at point r
     *  \param[out] de the value of the derivative of the dielectric constant
     *                 at point r
     *  \param[in]   r evaluation point
     */
    void operator()(double & e, double & de, const double r) const
    {
        e = value(r);
        de = derivative(r);
    }
    double epsilon1() const { return epsilon1_; }
    double epsilon2() const { return epsilon2_; }
    double width() const { return width_; }
    double center() const { return center_; }
    friend std::ostream & operator<<(std::ostream & os, TanhDiffuse & th) {
        return th.printObject(os);
    }
};

#endif // TANHDIFFUSE_HPP
