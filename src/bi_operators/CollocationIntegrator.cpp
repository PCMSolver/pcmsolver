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

#include "CollocationIntegrator.hpp"

#include <cmath>
#include <stdexcept>

#include "Config.hpp"

#include <Eigen/Dense>

#include "Element.hpp"

double CollocationIntegrator::computeS(const Vacuum<double> * /* gf */, const Element & e) const {
	double area = e.area();
	return (factor_ * std::sqrt(4 * M_PI / area));
}
double CollocationIntegrator::computeS(const Vacuum<AD_directional> * /* gf */, const Element & e) const {
	double area = e.area();
	return (factor_ * std::sqrt(4 * M_PI / area));
}
double CollocationIntegrator::computeS(const Vacuum<AD_gradient> * /* gf */, const Element & e) const {
	double area = e.area();
	return (factor_ * std::sqrt(4 * M_PI / area));
}
double CollocationIntegrator::computeS(const Vacuum<AD_hessian> * /* gf */, const Element & e) const {
	double area = e.area();
	return (factor_ * std::sqrt(4 * M_PI / area));
}

double CollocationIntegrator::computeD(const Vacuum<double> * /* gf */, const Element & e) const {
	double area = e.area();
	double radius = e.sphere().radius();
        return (-factor_ * std::sqrt(M_PI/ area) * (1.0 / radius));
}
double CollocationIntegrator::computeD(const Vacuum<AD_directional> * /* gf */, const Element & e) const {
	double area = e.area();
	double radius = e.sphere().radius();
        return (-factor_ * std::sqrt(M_PI/ area) * (1.0 / radius));
}
double CollocationIntegrator::computeD(const Vacuum<AD_gradient> * /* gf */, const Element & e) const {
	double area = e.area();
	double radius = e.sphere().radius();
        return (-factor_ * std::sqrt(M_PI/ area) * (1.0 / radius));
}
double CollocationIntegrator::computeD(const Vacuum<AD_hessian> * /* gf */, const Element & e) const {
	double area = e.area();
	double radius = e.sphere().radius();
        return (-factor_ * std::sqrt(M_PI/ area) * (1.0 / radius));
}

double CollocationIntegrator::computeS(const UniformDielectric<double> * gf, const Element & e) const {
	double epsInv = 1.0 / gf->epsilon();
	double area = e.area();
	return (factor_ * std::sqrt(4 * M_PI / area) * epsInv);
}
double CollocationIntegrator::computeS(const UniformDielectric<AD_directional> * gf, const Element & e) const {
	double epsInv = 1.0 / gf->epsilon();
	double area = e.area();
	return (factor_ * std::sqrt(4 * M_PI / area) * epsInv);
}
double CollocationIntegrator::computeS(const UniformDielectric<AD_gradient> * gf, const Element & e) const {
	double epsInv = 1.0 / gf->epsilon();
	double area = e.area();
	return (factor_ * std::sqrt(4 * M_PI / area) * epsInv);
}
double CollocationIntegrator::computeS(const UniformDielectric<AD_hessian> * gf, const Element & e) const {
	double epsInv = 1.0 / gf->epsilon();
	double area = e.area();
	return (factor_ * std::sqrt(4 * M_PI / area) * epsInv);
}

double CollocationIntegrator::computeD(const UniformDielectric<double> * /* gf */, const Element & e) const {
	double area = e.area();
	double radius = e.sphere().radius();
    return (-factor_ * std::sqrt(M_PI/ area) * (1.0 / radius));
}
double CollocationIntegrator::computeD(const UniformDielectric<AD_directional> * /* gf */, const Element & e) const {
	double area = e.area();
	double radius = e.sphere().radius();
    return (-factor_ * std::sqrt(M_PI/ area) * (1.0 / radius));
}
double CollocationIntegrator::computeD(const UniformDielectric<AD_gradient> * /* gf */, const Element & e) const {
	double area = e.area();
	double radius = e.sphere().radius();
    return (-factor_ * std::sqrt(M_PI/ area) * (1.0 / radius));
}
double CollocationIntegrator::computeD(const UniformDielectric<AD_hessian> * /* gf */, const Element & e) const {
	double area = e.area();
	double radius = e.sphere().radius();
    return (-factor_ * std::sqrt(M_PI/ area) * (1.0 / radius));
}

double CollocationIntegrator::computeS(const IonicLiquid<double> * /* gf */, const Element & /* e */) const {
	return 0.0;
}
double CollocationIntegrator::computeS(const IonicLiquid<AD_directional> * /* gf */, const Element & /* e */) const {
	return 0.0;
}
double CollocationIntegrator::computeS(const IonicLiquid<AD_gradient> * /* gf */, const Element & /* e */) const {
	return 0.0;
}
double CollocationIntegrator::computeS(const IonicLiquid<AD_hessian> * /* gf */, const Element & /* e */) const {
	return 0.0;
}

double CollocationIntegrator::computeD(const IonicLiquid<double> * /* gf */, const Element & /* e */) const {
	return 0.0;
}
double CollocationIntegrator::computeD(const IonicLiquid<AD_directional> * /* gf */, const Element & /* e */) const {
	return 0.0;
}
double CollocationIntegrator::computeD(const IonicLiquid<AD_gradient> * /* gf */, const Element & /* e */) const {
	return 0.0;
}
double CollocationIntegrator::computeD(const IonicLiquid<AD_hessian> * /* gf */, const Element & /* e */) const {
	return 0.0;
}

double CollocationIntegrator::computeS(const AnisotropicLiquid<double> * /* gf */, const Element & /* e */) const {
	return 0.0;
}
double CollocationIntegrator::computeS(const AnisotropicLiquid<AD_directional> * /* gf */, const Element & /* e */) const {
	return 0.0;
}
double CollocationIntegrator::computeS(const AnisotropicLiquid<AD_gradient> * /* gf */, const Element & /* e */) const {
	return 0.0;
}
double CollocationIntegrator::computeS(const AnisotropicLiquid<AD_hessian> * /* gf */, const Element & /* e */) const {
	return 0.0;
}

double CollocationIntegrator::computeD(const AnisotropicLiquid<double> * /* gf */, const Element & /* e */) const {
	return 0.0;
}
double CollocationIntegrator::computeD(const AnisotropicLiquid<AD_directional> * /* gf */, const Element & /* e */) const {
	return 0.0;
}
double CollocationIntegrator::computeD(const AnisotropicLiquid<AD_gradient> * /* gf */, const Element & /* e */) const {
	return 0.0;
}
double CollocationIntegrator::computeD(const AnisotropicLiquid<AD_hessian> * /* gf */, const Element & /* e */) const {
	return 0.0;
}

double CollocationIntegrator::computeS(const TanhSphericalDiffuse * gf, const Element & e) const {
    // The singular part is "integrated" as usual, while the nonsingular part is evaluated in full
	double area = e.area();
    // Diagonal of S inside the cavity
    double Sii_I = factor_ * std::sqrt(4 * M_PI / area);
    // "Diagonal" of Coulomb singularity separation coefficient
    double coulomb_coeff = gf->coefficientCoulomb(e.center(), e.center());
    // "Diagonal" of the image Green's function
    double image = gf->imagePotential(e.center(), e.center());

	return (Sii_I / coulomb_coeff + image);
}

double CollocationIntegrator::computeD(const TanhSphericalDiffuse * gf, const Element & e) const {
    // The singular part is "integrated" as usual, while the nonsingular part is evaluated in full
	double area = e.area();
	double radius = e.sphere().radius();
    // Diagonal of S inside the cavity
    double Sii_I = factor_ * std::sqrt(4 * M_PI / area);
    // Diagonal of D inside the cavity
    double Dii_I = -factor_ * std::sqrt(M_PI/ area) * (1.0 / radius);
    // "Diagonal" of Coulomb singularity separation coefficient
    double coulomb_coeff = gf->coefficientCoulomb(e.center(), e.center());
    // "Diagonal" of the directional derivative of the Coulomb singularity separation coefficient
    double coeff_grad = gf->coefficientCoulombDerivative(e.normal(), e.center(), e.center()) / std::pow(coulomb_coeff, 2);
    // "Diagonal" of the directional derivative of the image Green's function
    double image_grad = gf->imagePotentialDerivative(e.normal(), e.center(), e.center());

    double eps_r2 = 0.0;
    double d_eps_r2 = 0.0;
    gf->epsilon(eps_r2, d_eps_r2, e.center());

    return eps_r2 * (Dii_I / coulomb_coeff - Sii_I * coeff_grad + image_grad);
}
