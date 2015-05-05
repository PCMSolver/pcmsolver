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

#include "NumericalIntegrator.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>

#include "Config.hpp"

#include <Eigen/Dense>

#include "Element.hpp"
#include "MathUtils.hpp"
#include "QuadratureRules.hpp"
#include "Sphere.hpp"

using namespace std::placeholders;

typedef std::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)>
singleLayerIntegrand;

typedef std::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &,
		        const Eigen::Vector3d &)> doubleLayerIntegrand;

double NumericalIntegrator::computeS(const Vacuum<double> * gf, const Element & e) const {
	singleLayerIntegrand F = std::bind(&Vacuum<double>::function, *gf, _1, _2);
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeS(const Vacuum<AD_directional> * gf, const Element & e) const {
	singleLayerIntegrand F = std::bind(&Vacuum<AD_directional>::function, *gf, _1, _2);
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeS(const Vacuum<AD_gradient> * gf, const Element & e) const {
	singleLayerIntegrand F = std::bind(&Vacuum<AD_gradient>::function, *gf, _1, _2);
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeS(const Vacuum<AD_hessian> * gf, const Element & e) const {
	singleLayerIntegrand F = std::bind(&Vacuum<AD_hessian>::function, *gf, _1, _2);
        return integrator<32, 16>(F, e);
}

double NumericalIntegrator::computeD(const Vacuum<double> * gf, const Element & e) const {
	doubleLayerIntegrand F = std::bind(&Vacuum<double>::derivative, *gf, _1, _2, _3);
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeD(const Vacuum<AD_directional> * gf, const Element & e) const {
	doubleLayerIntegrand F = std::bind(&Vacuum<AD_directional>::derivative, *gf, _1, _2, _3);
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeD(const Vacuum<AD_gradient> * gf, const Element & e) const {
	doubleLayerIntegrand F = std::bind(&Vacuum<AD_gradient>::derivative, *gf, _1, _2, _3);
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeD(const Vacuum<AD_hessian> * gf, const Element & e) const {
	doubleLayerIntegrand F = std::bind(&Vacuum<AD_hessian>::derivative, *gf, _1, _2, _3);
        return integrator<32, 16>(F, e);
}

double NumericalIntegrator::computeS(const UniformDielectric<double> * gf, const Element & e) const {
	singleLayerIntegrand F = std::bind(&UniformDielectric<double>::function, *gf, _1, _2);
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeS(const UniformDielectric<AD_directional> * gf, const Element & e) const {
	singleLayerIntegrand F = std::bind(&UniformDielectric<AD_directional>::function, *gf, _1, _2);
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeS(const UniformDielectric<AD_gradient> * gf, const Element & e) const {
	singleLayerIntegrand F = std::bind(&UniformDielectric<AD_gradient>::function, *gf, _1, _2);
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeS(const UniformDielectric<AD_hessian> * gf, const Element & e) const {
	singleLayerIntegrand F = std::bind(&UniformDielectric<AD_hessian>::function, *gf, _1, _2);
        return integrator<32, 16>(F, e);
}

double NumericalIntegrator::computeD(const UniformDielectric<double> * gf, const Element & e) const {
	doubleLayerIntegrand F = std::bind(&UniformDielectric<double>::derivative, *gf, _1, _2, _3);
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeD(const UniformDielectric<AD_directional> * gf, const Element & e) const {
	doubleLayerIntegrand F = std::bind(&UniformDielectric<AD_directional>::derivative, *gf, _1, _2, _3);
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeD(const UniformDielectric<AD_gradient> * gf, const Element & e) const {
	doubleLayerIntegrand F = std::bind(&UniformDielectric<AD_gradient>::derivative, *gf, _1, _2, _3);
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeD(const UniformDielectric<AD_hessian> * gf, const Element & e) const {
	doubleLayerIntegrand F = std::bind(&UniformDielectric<AD_hessian>::derivative, *gf, _1, _2, _3);
        return integrator<32, 16>(F, e);
}

double NumericalIntegrator::computeS(const IonicLiquid<double> * gf, const Element & e) const {
	singleLayerIntegrand F = std::bind(&IonicLiquid<double>::function, *gf, _1, _2);
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeS(const IonicLiquid<AD_directional> * gf, const Element & e) const {
	singleLayerIntegrand F = std::bind(&IonicLiquid<AD_directional>::function, *gf, _1, _2);
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeS(const IonicLiquid<AD_gradient> * gf, const Element & e) const {
	singleLayerIntegrand F = std::bind(&IonicLiquid<AD_gradient>::function, *gf, _1, _2);
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeS(const IonicLiquid<AD_hessian> * gf, const Element & e) const {
	singleLayerIntegrand F = std::bind(&IonicLiquid<AD_hessian>::function, *gf, _1, _2);
        return integrator<32, 16>(F, e);
}

double NumericalIntegrator::computeD(const IonicLiquid<double> * gf, const Element & e) const {
	doubleLayerIntegrand F = std::bind(&IonicLiquid<double>::derivative, *gf, _1, _2, _3);
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeD(const IonicLiquid<AD_directional> * gf, const Element & e) const {
	doubleLayerIntegrand F = std::bind(&IonicLiquid<AD_directional>::derivative, *gf, _1, _2, _3);
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeD(const IonicLiquid<AD_gradient> * gf, const Element & e) const {
	doubleLayerIntegrand F = std::bind(&IonicLiquid<AD_gradient>::derivative, *gf, _1, _2, _3);
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeD(const IonicLiquid<AD_hessian> * gf, const Element & e) const {
	doubleLayerIntegrand F = std::bind(&IonicLiquid<AD_hessian>::derivative, *gf, _1, _2, _3);
        return integrator<32, 16>(F, e);
}

double NumericalIntegrator::computeS(const AnisotropicLiquid<double> * gf, const Element & e) const {
	singleLayerIntegrand F = std::bind(&AnisotropicLiquid<double>::function, *gf, _1, _2);
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeS(const AnisotropicLiquid<AD_directional> * gf, const Element & e) const {
	singleLayerIntegrand F = std::bind(&AnisotropicLiquid<AD_directional>::function, *gf, _1, _2);
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeS(const AnisotropicLiquid<AD_gradient> * gf, const Element & e) const {
	singleLayerIntegrand F = std::bind(&AnisotropicLiquid<AD_gradient>::function, *gf, _1, _2);
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeS(const AnisotropicLiquid<AD_hessian> * gf, const Element & e) const {
	singleLayerIntegrand F = std::bind(&AnisotropicLiquid<AD_hessian>::function, *gf, _1, _2);
        return integrator<32, 16>(F, e);
}

double NumericalIntegrator::computeD(const AnisotropicLiquid<double> * gf, const Element & e) const {
	doubleLayerIntegrand F = std::bind(&AnisotropicLiquid<double>::derivative, *gf, _1, _2, _3);
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeD(const AnisotropicLiquid<AD_directional> * gf, const Element & e) const {
	doubleLayerIntegrand F = std::bind(&AnisotropicLiquid<AD_directional>::derivative, *gf, _1, _2, _3);
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeD(const AnisotropicLiquid<AD_gradient> * gf, const Element & e) const {
	doubleLayerIntegrand F = std::bind(&AnisotropicLiquid<AD_gradient>::derivative, *gf, _1, _2, _3);
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeD(const AnisotropicLiquid<AD_hessian> * gf, const Element & e) const {
	doubleLayerIntegrand F = std::bind(&AnisotropicLiquid<AD_hessian>::derivative, *gf, _1, _2, _3);
        return integrator<32, 16>(F, e);
}

double NumericalIntegrator::computeS(const TanhSphericalDiffuse * /* gf */, const Element & /* e */) const {
    // The singular part is "integrated" as usual, while the nonsingular part is evaluated in full
	return 0.0;
}

double NumericalIntegrator::computeD(const TanhSphericalDiffuse * /* gf */, const Element & /* e */) const {
    // The singular part is "integrated" as usual, while the nonsingular part is evaluated in full
	return 0.0;
}
