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

double CollocationIntegrator::computeS(const Vacuum<double> * gf, double area) const {
	return (factor_ * std::sqrt(4 * M_PI / area));
} 
double CollocationIntegrator::computeS(const Vacuum<AD_directional> * gf, double area) const {
	return (factor_ * std::sqrt(4 * M_PI / area));
}
double CollocationIntegrator::computeS(const Vacuum<AD_gradient> * gf, double area) const {
	return (factor_ * std::sqrt(4 * M_PI / area));
}
double CollocationIntegrator::computeS(const Vacuum<AD_hessian> * gf, double area) const {
	return (factor_ * std::sqrt(4 * M_PI / area));
}

double CollocationIntegrator::computeD(const Vacuum<double> * gf, double area, double radius) const {
        return (-factor_ * std::sqrt(M_PI/ area) * (1.0 / radius));
}
double CollocationIntegrator::computeD(const Vacuum<AD_directional> * gf, double area, double radius) const {
        return (-factor_ * std::sqrt(M_PI/ area) * (1.0 / radius));
}
double CollocationIntegrator::computeD(const Vacuum<AD_gradient> * gf, double area, double radius) const {
        return (-factor_ * std::sqrt(M_PI/ area) * (1.0 / radius));
}
double CollocationIntegrator::computeD(const Vacuum<AD_hessian> * gf, double area, double radius) const {
        return (-factor_ * std::sqrt(M_PI/ area) * (1.0 / radius));
}

double CollocationIntegrator::computeS(const UniformDielectric<double> * gf, double area) const {
	double epsInv = 1.0 / gf->epsilon();
	return (factor_ * std::sqrt(4 * M_PI / area) * epsInv);
}
double CollocationIntegrator::computeS(const UniformDielectric<AD_directional> * gf, double area) const {
	double epsInv = 1.0 / gf->epsilon();
	return (factor_ * std::sqrt(4 * M_PI / area) * epsInv);
}
double CollocationIntegrator::computeS(const UniformDielectric<AD_gradient> * gf, double area) const {
	double epsInv = 1.0 / gf->epsilon();
	return (factor_ * std::sqrt(4 * M_PI / area) * epsInv);
}
double CollocationIntegrator::computeS(const UniformDielectric<AD_hessian> * gf, double area) const {
	double epsInv = 1.0 / gf->epsilon();
	return (factor_ * std::sqrt(4 * M_PI / area) * epsInv);
}

double CollocationIntegrator::computeD(const UniformDielectric<double> * gf, double area, double radius) const {
	double epsInv = 1.0 / gf->epsilon();
        return (-factor_ * std::sqrt(M_PI/ area) * (1.0 / radius) * epsInv);
}
double CollocationIntegrator::computeD(const UniformDielectric<AD_directional> * gf, double area, double radius) const {
	double epsInv = 1.0 / gf->epsilon();
        return (-factor_ * std::sqrt(M_PI/ area) * (1.0 / radius) * epsInv);
}
double CollocationIntegrator::computeD(const UniformDielectric<AD_gradient> * gf, double area, double radius) const {
	double epsInv = 1.0 / gf->epsilon();
        return (-factor_ * std::sqrt(M_PI/ area) * (1.0 / radius) * epsInv);
}
double CollocationIntegrator::computeD(const UniformDielectric<AD_hessian> * gf, double area, double radius) const {
	double epsInv = 1.0 / gf->epsilon();
        return (-factor_ * std::sqrt(M_PI/ area) * (1.0 / radius) * epsInv);
}

double CollocationIntegrator::computeS(const IonicLiquid<double> * gf, double area) const {
	return 0.0;
}
double CollocationIntegrator::computeS(const IonicLiquid<AD_directional> * gf, double area) const {
	return 0.0;
}
double CollocationIntegrator::computeS(const IonicLiquid<AD_gradient> * gf, double area) const {
	return 0.0;
}
double CollocationIntegrator::computeS(const IonicLiquid<AD_hessian> * gf, double area) const {
	return 0.0;
}

double CollocationIntegrator::computeD(const IonicLiquid<double> * gf, double area, double radius) const {
	return 0.0;
}
double CollocationIntegrator::computeD(const IonicLiquid<AD_directional> * gf, double area, double radius) const {
	return 0.0;
}
double CollocationIntegrator::computeD(const IonicLiquid<AD_gradient> * gf, double area, double radius) const {
	return 0.0;
}
double CollocationIntegrator::computeD(const IonicLiquid<AD_hessian> * gf, double area, double radius) const {
	return 0.0;
}
