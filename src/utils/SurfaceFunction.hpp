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

#ifndef SURFACEFUNCTION_HPP
#define SURFACEFUNCTION_HPP

#include <iosfwd>
#include <string>

#include "Config.hpp"

#include <Eigen/Dense>

/*!
 * \file SurfaceFunction.hpp
 * \class SurfaceFunction
 * \brief A basic surface function class
 * \author Luca Frediani, Roberto Di Remigio
 * \date 2012, 2013
 *
 * This class is basically a wrapper around vectors containing electrostatic potentials
 * and apparent charges. Use judiciously, i.e. DO NOT use it directly in the core
 * classes (cavities, solvers) to avoid high coupling.
 */


class SurfaceFunction
{
public:
    SurfaceFunction() : name_(""), nPoints_(0), allocated_(false) {}
    SurfaceFunction(const std::string & n) 
	    : name_(n), nPoints_(0), allocated_(false) {}
    SurfaceFunction(const std::string & n, int np) 
	    : name_(n), nPoints_(np) {
        values_ = Eigen::VectorXd::Zero(nPoints_);
        allocated_ = true;
    }
    SurfaceFunction(const std::string & n, int np, double * v) 
	    : name_(n), nPoints_(np) {
        values_ = Eigen::VectorXd::Zero(nPoints_);
        allocated_ = true;
        Eigen::Map<Eigen::VectorXd> mapped_data(v, nPoints_);
	values_ = mapped_data;
    }
    SurfaceFunction(const std::string & n, int np, const Eigen::VectorXd & v) 
	    : name_(n), nPoints_(np), values_(v), allocated_(true) {} 
    ~SurfaceFunction() {
        allocated_ = false;
    }

    /// Copy constructor
    SurfaceFunction(const SurfaceFunction & other) 
	    : name_(other.name_), nPoints_(other.nPoints_), values_(other.values_) {
        allocated_ = true;
    }

    friend inline void swap(SurfaceFunction & left, SurfaceFunction & right);
    inline void swap(SurfaceFunction & other);
    /// Assignment operator.
    SurfaceFunction & operator=(SurfaceFunction other);
    /// Multiplication operator: product of two SurfaceFunctions version (scalar product of the values vectors).
    double operator*(const SurfaceFunction & other) const;
    /// Addition-assignment operator.
    SurfaceFunction & operator+=(const SurfaceFunction & other);
    /// Subtraction-assignment operator.
    SurfaceFunction & operator-=(const SurfaceFunction & other);
    /// Multiplication-assignment operator. Defined only for the uniform scaling case.
    SurfaceFunction & operator*=(double scaling);
    /// Division-assignment operator. Defined only for the uniform scaling case.
    SurfaceFunction & operator/=(double scaling);

    std::string & name() { return name_; }
    int nPoints() const { return nPoints_; }
    void value(int index, double value) { values_(index) = value; }
    double value(int index) const { return values_(index); }
    Eigen::VectorXd & vector() { return values_; }
    const Eigen::VectorXd & vector() const { return values_; }
    void allocate(int np) { values_.resize(np); }
    bool allocated() const { return allocated_; }
    void clear();

    void setValues(double * v);
    void getValues(double * v);

    friend std::ostream & operator<<(std::ostream & os, SurfaceFunction & sf) {
        return sf.printObject(os);
    }

private:
    std::ostream & printObject(std::ostream & os);
    std::string name_;
    int nPoints_;
    Eigen::VectorXd values_;
    bool allocated_;
};

/*!
 * \fn inline SurfaceFunction operator+(SurfaceFunction left, const SurfaceFunction & right)
 * \brief Addition operator
 * \param left the left hand side of the addition
 * \param right the right hand side of the addition
 */
inline SurfaceFunction operator+(SurfaceFunction left, const SurfaceFunction & right)
{
    left += right;
    return left;
}

/*!
 * \fn inline SurfaceFunction operator-(SurfaceFunction left, const SurfaceFunction & right)
 * \brief Subtraction operator
 * \param left the left hand side of the subtraction
 * \param right the right hand side of the subtraction
 */
inline SurfaceFunction operator-(SurfaceFunction left, const SurfaceFunction & right)
{
    left -= right;
    return left;
}

/*!
 * \fn inline SurfaceFunction operator*(double scaling, const SurfaceFunction & right)
 * \brief Multiplication operator: uniform scaling of SurfaceFunction version
 * \param scaling the scaling factor
 * \param object the surface function to be scaled
 */
inline SurfaceFunction operator*(double scaling, SurfaceFunction & object)
{
    object *= scaling;
    return object;
}

#endif // SURFACEFUNCTION_HPP
