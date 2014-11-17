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

#include "SurfaceFunction.hpp"

#include <string>
#include <stdexcept>

#include "Config.hpp"

#include <Eigen/Dense>

#include <boost/lexical_cast.hpp>

inline void swap(SurfaceFunction & left, SurfaceFunction & right)
{
    using std::swap;
    swap(left.name_, right.name_);
    swap(left.nPoints_,
         right.nPoints_); // This is maybe redundant as nPoints must be the same...
    swap(left.values_, right.values_);
}

inline void SurfaceFunction::swap(SurfaceFunction & other)
{
    using std::swap;
    swap(this->name_, other.name_);
    swap(this->nPoints_, other.nPoints_); // This is maybe redundant
    swap(this->values_, other.values_);
}

SurfaceFunction & SurfaceFunction::operator=(SurfaceFunction other)
{
    if (this == &other) // Check for self-assignment
        throw std::runtime_error("Are you trying to self-assign?");
    // No check on dimensions! This is assignment not comparison!!

    ::swap(*this, other);
    return *this;
}

double SurfaceFunction::operator*(const SurfaceFunction & other)
{
    if (this->nPoints_ != other.nPoints_)
        throw std::runtime_error("Incoherent dimensions of left and right operands!");
    return this->values_.dot(other.values_);
}

SurfaceFunction & SurfaceFunction::operator+=(const SurfaceFunction & other)
{
    if (this->nPoints_ != other.nPoints_)
        throw std::runtime_error("Incoherent dimensions of left and right operands!");
    this->name_ += "+" + other.name_;
    this->values_ = this->values_ + other.values_;
    return *this;
}

SurfaceFunction & SurfaceFunction::operator-=(const SurfaceFunction & other)
{
    if (this->nPoints_ != other.nPoints_)
        throw std::runtime_error("Incoherent dimensions of left and right operands!");
    this->name_ += "-" + other.name_;
    this->values_ = this->values_ - other.values_;
    return *this;
}

SurfaceFunction & SurfaceFunction::operator*=(double scaling)
{
    this->name_ = boost::lexical_cast<std::string>(scaling) + "*" + this->name_;
    this->values_ *= scaling;
    return *this;
}

SurfaceFunction & SurfaceFunction::operator/=(double scaling)
{
    if (scaling == 0.0)
        throw std::runtime_error("You are dividing by zero!");
    this->name_ = boost::lexical_cast<std::string>(scaling) + "/" + this->name_;
    this->values_ /= scaling;
    return *this;
}

void SurfaceFunction::setValues(double * v)
{
    if (!allocated_)
        throw std::runtime_error("Surface function not allocated!");
    // Zero out any previous value
    values_.setZero();

    for (int i = 0; i < nPoints_; ++i) {
        values_(i) = v[i];
    }
}

void SurfaceFunction::getValues(double * v)
{
    for (int i = 0; i < nPoints_; ++i) {
        v[i] = values_(i);
    }
}

void SurfaceFunction::clear()
{
    values_.setZero();
}

std::ostream & SurfaceFunction::printObject(std::ostream & os)
{
    os << "Surface Function " << name_ << std::endl;
    if (!allocated_)
        throw std::runtime_error("Surface function not allocated!");

    os << values_.transpose();

    return os;
}
