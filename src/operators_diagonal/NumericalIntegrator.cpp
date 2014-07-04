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

#include "NumericalIntegrator.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include <boost/bind.hpp>
#include <boost/bind/placeholders.hpp>
#include <boost/function.hpp>

#include "Element.hpp"
#include "MathUtils.hpp"
#include "QuadratureRules.hpp"
#include "Sphere.hpp"

typedef boost::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)>
singleLayerIntegrand;

typedef boost::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &, 
		        const Eigen::Vector3d &)> doubleLayerIntegrand;

double NumericalIntegrator::computeS(const Vacuum<double> * gf, const Element & e) const {
	singleLayerIntegrand F = boost::bind(&Vacuum<double>::function, *gf, _1, _2); 
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeS(const Vacuum<AD_directional> * gf, const Element & e) const {
	singleLayerIntegrand F = boost::bind(&Vacuum<AD_directional>::function, *gf, _1, _2); 
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeS(const Vacuum<AD_gradient> * gf, const Element & e) const {
	singleLayerIntegrand F = boost::bind(&Vacuum<AD_gradient>::function, *gf, _1, _2); 
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeS(const Vacuum<AD_hessian> * gf, const Element & e) const {
	singleLayerIntegrand F = boost::bind(&Vacuum<AD_hessian>::function, *gf, _1, _2); 
        return integrator<32, 16>(F, e);
}

double NumericalIntegrator::computeD(const Vacuum<double> * gf, const Element & e) const {
	doubleLayerIntegrand F = boost::bind(&Vacuum<double>::derivativeProbe, *gf, _1, _2, _3); 
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeD(const Vacuum<AD_directional> * gf, const Element & e) const {
	doubleLayerIntegrand F = boost::bind(&Vacuum<AD_directional>::derivativeProbe, *gf, _1, _2, _3); 
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeD(const Vacuum<AD_gradient> * gf, const Element & e) const {
	doubleLayerIntegrand F = boost::bind(&Vacuum<AD_gradient>::derivativeProbe, *gf, _1, _2, _3); 
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeD(const Vacuum<AD_hessian> * gf, const Element & e) const {
	doubleLayerIntegrand F = boost::bind(&Vacuum<AD_hessian>::derivativeProbe, *gf, _1, _2, _3); 
        return integrator<32, 16>(F, e);
}

double NumericalIntegrator::computeS(const UniformDielectric<double> * gf, const Element & e) const {
	singleLayerIntegrand F = boost::bind(&UniformDielectric<double>::function, *gf, _1, _2); 
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeS(const UniformDielectric<AD_directional> * gf, const Element & e) const {
	singleLayerIntegrand F = boost::bind(&UniformDielectric<AD_directional>::function, *gf, _1, _2); 
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeS(const UniformDielectric<AD_gradient> * gf, const Element & e) const {
	singleLayerIntegrand F = boost::bind(&UniformDielectric<AD_gradient>::function, *gf, _1, _2); 
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeS(const UniformDielectric<AD_hessian> * gf, const Element & e) const {
	singleLayerIntegrand F = boost::bind(&UniformDielectric<AD_hessian>::function, *gf, _1, _2); 
        return integrator<32, 16>(F, e);
}

double NumericalIntegrator::computeD(const UniformDielectric<double> * gf, const Element & e) const {
	doubleLayerIntegrand F = boost::bind(&UniformDielectric<double>::derivativeProbe, *gf, _1, _2, _3); 
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeD(const UniformDielectric<AD_directional> * gf, const Element & e) const {
	doubleLayerIntegrand F = boost::bind(&UniformDielectric<AD_directional>::derivativeProbe, *gf, _1, _2, _3); 
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeD(const UniformDielectric<AD_gradient> * gf, const Element & e) const {
	doubleLayerIntegrand F = boost::bind(&UniformDielectric<AD_gradient>::derivativeProbe, *gf, _1, _2, _3); 
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeD(const UniformDielectric<AD_hessian> * gf, const Element & e) const {
	doubleLayerIntegrand F = boost::bind(&UniformDielectric<AD_hessian>::derivativeProbe, *gf, _1, _2, _3); 
        return integrator<32, 16>(F, e);
}

double NumericalIntegrator::computeS(const IonicLiquid<double> * gf, const Element & e) const {
	singleLayerIntegrand F = boost::bind(&IonicLiquid<double>::function, *gf, _1, _2); 
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeS(const IonicLiquid<AD_directional> * gf, const Element & e) const {
	singleLayerIntegrand F = boost::bind(&IonicLiquid<AD_directional>::function, *gf, _1, _2); 
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeS(const IonicLiquid<AD_gradient> * gf, const Element & e) const {
	singleLayerIntegrand F = boost::bind(&IonicLiquid<AD_gradient>::function, *gf, _1, _2); 
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeS(const IonicLiquid<AD_hessian> * gf, const Element & e) const {
	singleLayerIntegrand F = boost::bind(&IonicLiquid<AD_hessian>::function, *gf, _1, _2); 
        return integrator<32, 16>(F, e);
}

double NumericalIntegrator::computeD(const IonicLiquid<double> * gf, const Element & e) const {
	doubleLayerIntegrand F = boost::bind(&IonicLiquid<double>::derivativeProbe, *gf, _1, _2, _3); 
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeD(const IonicLiquid<AD_directional> * gf, const Element & e) const {
	doubleLayerIntegrand F = boost::bind(&IonicLiquid<AD_directional>::derivativeProbe, *gf, _1, _2, _3); 
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeD(const IonicLiquid<AD_gradient> * gf, const Element & e) const {
	doubleLayerIntegrand F = boost::bind(&IonicLiquid<AD_gradient>::derivativeProbe, *gf, _1, _2, _3); 
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeD(const IonicLiquid<AD_hessian> * gf, const Element & e) const {
	doubleLayerIntegrand F = boost::bind(&IonicLiquid<AD_hessian>::derivativeProbe, *gf, _1, _2, _3); 
        return integrator<32, 16>(F, e);
}

double NumericalIntegrator::computeS(const AnisotropicLiquid<double> * gf, const Element & e) const {
	singleLayerIntegrand F = boost::bind(&AnisotropicLiquid<double>::function, *gf, _1, _2); 
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeS(const AnisotropicLiquid<AD_directional> * gf, const Element & e) const {
	singleLayerIntegrand F = boost::bind(&AnisotropicLiquid<AD_directional>::function, *gf, _1, _2); 
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeS(const AnisotropicLiquid<AD_gradient> * gf, const Element & e) const {
	singleLayerIntegrand F = boost::bind(&AnisotropicLiquid<AD_gradient>::function, *gf, _1, _2); 
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeS(const AnisotropicLiquid<AD_hessian> * gf, const Element & e) const {
	singleLayerIntegrand F = boost::bind(&AnisotropicLiquid<AD_hessian>::function, *gf, _1, _2); 
        return integrator<32, 16>(F, e);
}

double NumericalIntegrator::computeD(const AnisotropicLiquid<double> * gf, const Element & e) const {
	doubleLayerIntegrand F = boost::bind(&AnisotropicLiquid<double>::derivativeProbe, *gf, _1, _2, _3); 
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeD(const AnisotropicLiquid<AD_directional> * gf, const Element & e) const {
	doubleLayerIntegrand F = boost::bind(&AnisotropicLiquid<AD_directional>::derivativeProbe, *gf, _1, _2, _3); 
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeD(const AnisotropicLiquid<AD_gradient> * gf, const Element & e) const {
	doubleLayerIntegrand F = boost::bind(&AnisotropicLiquid<AD_gradient>::derivativeProbe, *gf, _1, _2, _3); 
        return integrator<32, 16>(F, e);
}
double NumericalIntegrator::computeD(const AnisotropicLiquid<AD_hessian> * gf, const Element & e) const {
	doubleLayerIntegrand F = boost::bind(&AnisotropicLiquid<AD_hessian>::derivativeProbe, *gf, _1, _2, _3); 
        return integrator<32, 16>(F, e);
}
