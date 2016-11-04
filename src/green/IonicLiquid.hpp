/**
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2016 Roberto Di Remigio, Luca Frediani and collaborators.
 *
 * This file is part of PCMSolver.
 *
 * PCMSolver is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PCMSolver is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to the
 * PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
 */

#ifndef IONICLIQUID_HPP
#define IONICLIQUID_HPP

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Core>

class Element;
struct Yukawa;

#include "DerivativeTypes.hpp"
#include "GreensFunction.hpp"

/*! \file IonicLiquid.hpp
 *  \class IonicLiquid
 *  \brief Green's functions for ionic liquid, described by the linearized
 * Poisson-Boltzmann equation.
 *  \author Luca Frediani, Roberto Di Remigio
 *  \date 2013-2016
 *  \tparam DerivativeTraits evaluation strategy for the function and its derivatives
 */

template <typename DerivativeTraits = AD_directional>
class IonicLiquid __final : public GreensFunction<DerivativeTraits, Yukawa> {
public:
  IonicLiquid(double eps, double k);
  virtual ~IonicLiquid() {}

  friend std::ostream & operator<<(std::ostream & os, IonicLiquid & gf) {
    return gf.printObject(os);
  }

private:
  virtual DerivativeTraits operator()(DerivativeTraits * sp,
                                      DerivativeTraits * pp) const __override;
  virtual double kernelD_impl(const Eigen::Vector3d & direction,
                              const Eigen::Vector3d & p1,
                              const Eigen::Vector3d & p2) const __override;

  virtual KernelS exportKernelS_impl() const __override;
  virtual KernelD exportKernelD_impl() const __override;

  virtual double singleLayer_impl(const Element & /* e */,
                                  double /* factor */) const __override;
  virtual double doubleLayer_impl(const Element & /* e */,
                                  double /* factor */) const __override;

  virtual std::ostream & printObject(std::ostream & os) __override;
};

#endif // IONICLIQUID_HPP
