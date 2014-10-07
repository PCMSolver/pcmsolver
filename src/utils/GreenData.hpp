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

#ifndef GREENDATA_HPP
#define GREENDATA_HPP

#include <vector>

#include "Config.hpp"

class DiagonalIntegrator;

/*! @struct greenData
 *  @brief Contains all data defined from user input in the green section.
 *  @var greenData::how
 *  The way the derivatives of the Green's function are evaluated.
 *  @var greenData::epsilon
 *  The permittivity.
 *  @var greenData::kappa
 *  Inverse of the Debye length.
 *  @var greenData::epsilonTensor
 *  Diagonal values of the permittivity tensor with respect to the lab frame
 *  @var greenData::eulertAngles
 *  Euler angles giving the rotation of the solvent orientation with respect to the lab frame
 *  Default is zero degrees for all three angles
 *  @var greenData::epsilonReal
 *  Real part of the permittivity of a metal sphere.
 *  @var greenData::epsilonImaginary
 *  Imaginary part of the permittivity of a metal sphere.
 *  @var greenData::NPspheres
 *  Coordinates of the metal sphere center.
 *  @var greenData::NPradii
 *  Radius of the the metal sphere.
 *  @var greenData::integrator
 *  strategy to calculate the diagonal elements of S and D operator
 */

struct greenData {
    int how;
    double epsilon;
    double kappa;
    Eigen::Vector3d epsilonTensor;
    Eigen::Vector3d eulerAngles;
    double epsilonReal;
    double epsilonImaginary;
    std::vector<double> NPspheres;
    double NPradii;
    DiagonalIntegrator * integrator;

    greenData(int _how, double _epsilon = 1.0, DiagonalIntegrator * _integrator = NULL, double _kappa = 0.0,
	      const Eigen::Vector3d & epstens = Eigen::Vector3d::Zero(), const Eigen::Vector3d & euler = Eigen::Vector3d::Zero(),
              double _epsReal = 0.0, double _epsImaginary = 0.0,
              const std::vector<double> & _sphere = std::vector<double>(),
              double _sphRadius = 0.0) :
        how(_how), epsilon(_epsilon), kappa(_kappa), epsilonTensor(epstens), eulerAngles(euler),
	epsilonReal(_epsReal), epsilonImaginary(_epsImaginary),
        NPspheres(_sphere), NPradii(_sphRadius), integrator(_integrator) {}
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW /* See http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html */
};

#endif // GREENDATA_HPP
