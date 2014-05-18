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

//#include <boost/numeric/quadrature/kronrodgauss.hpp>
//#include <boost/numeric/quadrature/error_estimator.hpp>

#include "Element.hpp"
#include "MathUtils.hpp"
#include "Sphere.hpp"

double NumericalIntegrator::computeS(const Vacuum<double> * gf, const Element & e) const {
	return 0.0;
} 
double NumericalIntegrator::computeS(const Vacuum<AD_directional> * gf, const Element & e) const {
	return 0.0;
}
double NumericalIntegrator::computeS(const Vacuum<AD_gradient> * gf, const Element & e) const {
	return 0.0;
}
double NumericalIntegrator::computeS(const Vacuum<AD_hessian> * gf, const Element & e) const {
	return 0.0;
}

double NumericalIntegrator::computeD(const Vacuum<double> * gf, const Element & e) const {
	return 0.0;
}
double NumericalIntegrator::computeD(const Vacuum<AD_directional> * gf, const Element & e) const {
	return 0.0;
}
double NumericalIntegrator::computeD(const Vacuum<AD_gradient> * gf, const Element & e) const {
	return 0.0;
}
double NumericalIntegrator::computeD(const Vacuum<AD_hessian> * gf, const Element & e) const {
	return 0.0;
}

double NumericalIntegrator::computeS(const UniformDielectric<double> * gf, const Element & e) const {
	return 0.0;
}
double NumericalIntegrator::computeS(const UniformDielectric<AD_directional> * gf, const Element & e) const {
	return 0.0;
}
double NumericalIntegrator::computeS(const UniformDielectric<AD_gradient> * gf, const Element & e) const {
	return 0.0;
}
double NumericalIntegrator::computeS(const UniformDielectric<AD_hessian> * gf, const Element & e) const {
	return 0.0;
}

double NumericalIntegrator::computeD(const UniformDielectric<double> * gf, const Element & e) const {
	return 0.0;
}
double NumericalIntegrator::computeD(const UniformDielectric<AD_directional> * gf, const Element & e) const {
	return 0.0;
}
double NumericalIntegrator::computeD(const UniformDielectric<AD_gradient> * gf, const Element & e) const {
	return 0.0;
}
double NumericalIntegrator::computeD(const UniformDielectric<AD_hessian> * gf, const Element & e) const {
	return 0.0;
}

double NumericalIntegrator::computeS(const IonicLiquid<double> * gf, const Element & e) const {
	return 0.0;
}
double NumericalIntegrator::computeS(const IonicLiquid<AD_directional> * gf, const Element & e) const {
	return 0.0;
}
double NumericalIntegrator::computeS(const IonicLiquid<AD_gradient> * gf, const Element & e) const {
	return 0.0;
}
double NumericalIntegrator::computeS(const IonicLiquid<AD_hessian> * gf, const Element & e) const {
	return 0.0;
}

double NumericalIntegrator::computeD(const IonicLiquid<double> * gf, const Element & e) const {
	return 0.0;
}
double NumericalIntegrator::computeD(const IonicLiquid<AD_directional> * gf, const Element & e) const {
	return 0.0;
}
double NumericalIntegrator::computeD(const IonicLiquid<AD_gradient> * gf, const Element & e) const {
	return 0.0;
}
double NumericalIntegrator::computeD(const IonicLiquid<AD_hessian> * gf, const Element & e) const {
	return 0.0;
}

double NumericalIntegrator::computeS(const AnisotropicLiquid<double> * gf, const Element & e) const {
	double Sii = 0.0, S64 = 0.0;
	// Extract relevant data from Element
	int nVertices = e.nVertices();
	double area = e.area();
	Eigen::Vector3d center = e.center();
	Eigen::Vector3d normal = e.normal();
	Sphere sph = e.sphere();
	Eigen::Matrix3Xd vertices = e.vertices();
	Eigen::Matrix3Xd arcs = e.arcs();

	double xx = 0.0, 
	       xy = 0.0, 
	       xz = normal(0), 
	       yy = 0.0, 
	       yx = 0.0, 
	       yz = normal(1), 
	       zz = normal(2), 
	       zx = 0.0, 
	       zy = 0.0;

	double rmin = 0.99; // Some kind of threshold
	if (std::abs(xz) <= rmin) {
		rmin = std::abs(xz);
		xx = 0.0;
		yx = -zz / std::sqrt(1.0 - std::pow(xz, 2));
		zx =  yz / std::sqrt(1.0 - std::pow(xz, 2));
	}
	if (std::abs(yz) <= rmin) {
		rmin = std::abs(yz);
		xx =  yz / std::sqrt(1.0 - std::pow(yz, 2));
		yx = 0.0;
		zx = -xz / std::sqrt(1.0 - std::pow(yz, 2));
	}
	if (std::abs(zz) <= rmin) {
		rmin = std::abs(yz);
		xx =  yz / std::sqrt(1.0 - std::pow(zz, 2));
		yx = -xz / std::sqrt(1.0 - std::pow(zz, 2));
		zx = 0.0;
	}
	Eigen::Vector3d tangent; 
	tangent << xx, yx, zx; // Not sure this is the tangent vector...
	
	/*xy = yz * zx - yx * zz;
	yy = zz * xx - zx * xz;
	zy = xz * yx - xx * yz;
	Eigen::Vector3d bitangent << xy, yy, zy; // Not sure this is the bitangent vector...*/
	Eigen::Vector3d bitangent = normal.cross(tangent);

	std::vector<double> theta(nVertices), phi(nVertices), phinumb(nVertices+1);
	std::vector<int> numb(nVertices+1);
	// Clean-up heap crap
	std::fill_n(theta.begin(),   nVertices,   0.0);
	std::fill_n(phi.begin(),     nVertices,   0.0);
	std::fill_n(numb.begin(),    nVertices+1, 0);
	std::fill_n(phinumb.begin(), nVertices+1, 0.0);
        // Calculate a number of angles
	for (int i = 0; i < nVertices; ++i) {
		Eigen::Vector3d vertex_normal = (vertices.col(i) - sph.center()) / sph.radius();
		double scal1 = vertex_normal.dot(normal);
		if (scal1 >=  1.0) scal1 = 1.0;
		if (scal1 <= -1.0) scal1 = -1.0;
		theta[i] = std::acos(scal1);
		double scal2 = vertex_normal.dot(tangent) / std::sin(theta[i]);
		if (scal2 >=  1.0) scal2 = 1.0;
		if (scal2 <= -1.0) scal2 = -1.0;
		phi[i] = std::acos(scal2);
		double sin_phi = vertex_normal.dot(bitangent) / std::sin(theta[i]);
		if (sin_phi <= 0.0) phi[i] = 2 * M_PI - phi[i];
		if (i != 0) {
			phi[i] -= phi[0];
			if (phi[i] < 0.0) phi[i] += 2 * M_PI;
		}
	}
	// Recalculate tangent and bitangent vectors
	tangent *= std::cos(phi[0]);
	tangent += bitangent * std::sin(phi[0]);
	bitangent = normal.cross(tangent);
       
	// Populate numb and phinumb arrays
	phi[0] = 0.0;
	numb[0] = 1; numb[1] = 2;
	phinumb[0] = phi[0]; phinumb[1] = phinumb[1];
	for (int i = 2; i < nVertices; ++i) { // This loop is 2-based
		for (int j = 1; j < (i - 1); ++j) { // This loop is 1-based 
			if (phi[i] < phinumb[j]) {
				for (int k = 0; k < (i - j); ++k) {
					numb[i - k + 1] = numb[i - k];
					phinumb[i - k + 1] = phinumb[i - k];
				}
				numb[j] = i;
				phinumb[j] = phi[i];
				goto jump; // Ugly, to be refactored!!!
			}
			numb[i] = i;
			phinumb[i] = phi[i];
		}
		jump:
		std::cout << "jumped" << std::endl;
	}
	numb[nVertices] = numb[0];
	phinumb[nVertices] = 2 * M_PI;

	for (int i = 0; i < nVertices; ++i) { // Loop on edges
		double pha = phinumb[i];
		double phb = phinumb[i+1];
		double aph = (pha - phb) / 2.0;
		double bph = (pha + phb) / 2.0;
		double tha = theta[numb[i]];
		double thb = theta[numb[i+1]];
		double thmax = 0.0;
		Eigen::Vector3d oc = (arcs.col(i) - sph.center()) / sph.radius();
		double oc_norm = oc.norm();
		double oc_norm2 = std::pow(oc_norm, 2);
		for (int j = 0; j < 32; ++j) { // Loop on Gaussian points (64-points rule)
			for (int k = 0; k <= 1; ++k) {
				double ph = (2*k - 1) * aph * gauss64Abscissa(j) + bph;
				double cos_ph = std::cos(ph);
				double sin_ph = std::sin(ph);
				if (oc_norm2 < 1.0e-07) {
					double cotg_thmax = std::sin(ph-pha) / std::tan(thb) + std::sin(phb-ph) / std::tan(tha);
					cotg_thmax /= std::sin(phb-pha);
					thmax = std::atan(1.0 / cotg_thmax);
				} else {
					Eigen::Vector3d scratch; 
					scratch << tangent.dot(oc), bitangent.dot(oc), normal.dot(oc);
					double aa = std::pow(tangent.dot(oc)*cos_ph + bitangent.dot(oc)*sin_ph, 2) + std::pow(normal.dot(oc), 2);
					double bb = -normal.dot(oc) * oc_norm2;
					double cc = std::pow(oc_norm2, 2) - std::pow(tangent.dot(oc)*cos_ph + bitangent.dot(oc)*sin_ph, 2);
					double ds = std::pow(bb, 2) - aa*cc;
					if (ds < 0.0) ds = 0.0;
					double cs = (-bb + std::sqrt(ds)) / aa;
					if (cs > 1.0) cs = 1.0;
					if (cs < -1.0) cs = 1.0;
					thmax = std::acos(cs);
				}
				double ath = thmax / 2.0;

				double S16 = 0.0;
				if (thmax < 1.0e-08) goto jump1; // Ugly, to be refactored!!!
				for (int l = 0; l < 8; ++l) { // Loop on Gaussian points (16-points rule)
					for (int m = 0; m <= 1; ++m) {
						double th = (2*m - 1) * ath * gauss16Abscissa(l) + ath;
						double cos_th = std::cos(th);
						double sin_th = std::sin(th);
						double vx = tangent(0) * sin_th * cos_ph 
							  + bitangent(0) * sin_th * sin_ph 
							  + normal(0) * (cos_th - 1.0);
						double vy = tangent(1) * sin_th * cos_ph 
							  + bitangent(1) * sin_th * sin_ph 
							  + normal(1) * (cos_th - 1.0);
						double vz = tangent(2) * sin_th * cos_ph 
							  + bitangent(2) * sin_th * sin_ph 
							  + normal(2) * (cos_th - 1.0);
						Eigen::Vector3d evaluation; 
						evaluation << vx, vy, vz;
						double rth = std::sqrt(2 * (1.0 - cos_th));
						double rtheps = gf->function(evaluation, Eigen::Vector3d::Zero()); // Evaluate Green's function at Gaussian points
						S16 += sph.radius() * rtheps * sin_th * ath * gauss16Weight(l);
					}
				}
				S64 += S16 * aph * gauss64Weight(j);

				jump1:
				std::cout << "jumped!" << std::endl;
			}
		}
	}

	Sii = S64 * area;

        return Sii;
}
double NumericalIntegrator::computeS(const AnisotropicLiquid<AD_directional> * gf, const Element & e) const {
	return 0.0;
}
double NumericalIntegrator::computeS(const AnisotropicLiquid<AD_gradient> * gf, const Element & e) const {
	return 0.0;
}
double NumericalIntegrator::computeS(const AnisotropicLiquid<AD_hessian> * gf, const Element & e) const {
	return 0.0;
}

double NumericalIntegrator::computeD(const AnisotropicLiquid<double> * gf, const Element & e) const {
	return 0.0;
}
double NumericalIntegrator::computeD(const AnisotropicLiquid<AD_directional> * gf, const Element & e) const {
	return 0.0;
}
double NumericalIntegrator::computeD(const AnisotropicLiquid<AD_gradient> * gf, const Element & e) const {
	return 0.0;
}
double NumericalIntegrator::computeD(const AnisotropicLiquid<AD_hessian> * gf, const Element & e) const {
	return 0.0;
}
