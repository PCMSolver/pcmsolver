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
	double Sii = 0.0, S64 = 0.0;
	// Extract relevant data from Element
	int nVertices = e.nVertices();
	double area = e.area();
	Eigen::Vector3d center = e.center();
	Eigen::Vector3d normal = e.normal();
	Sphere sph = e.sphere();
	Eigen::Matrix3Xd vertices = e.vertices();
	Eigen::Matrix3Xd arcs = e.arcs();

	// Calculation of the tangent and the bitangent (binormal) vectors
	Eigen::Vector3d tangent, bitangent;
	e.tangent_and_bitangent(tangent, bitangent);

	std::vector<double> theta(nVertices), phi(nVertices), phinumb(nVertices+1);
	std::vector<int> numb(nVertices+1);
	// Clean-up heap crap
	std::fill_n(theta.begin(),   nVertices,   0.0);
	std::fill_n(phi.begin(),     nVertices,   0.0);
	std::fill_n(numb.begin(),    nVertices+1, 0);
	std::fill_n(phinumb.begin(), nVertices+1, 0.0);

        // Calculate the azimuthal and polar angles for the tessera vertices:
	// we use the normal, tangent and bitangent as a local reference frame
	for (int i = 0; i < nVertices; ++i) {
		Eigen::Vector3d vertex_normal = (vertices.col(i) - sph.center()) / sph.radius();
		// The cosine of the polar angle is given as the dot product of the normal at the vertex and the
		// normal at the tessera center: R\cos\theta
		double cos_theta = vertex_normal.dot(normal);
		if (cos_theta >=  1.0) cos_theta = 1.0;
		if (cos_theta <= -1.0) cos_theta = -1.0;
		theta[i] = std::acos(cos_theta);
		// The cosine of the azimuthal angle is given as the dot product of the normal at the vertex and the
		// tangent at the tessera center divided by the sine of the polar angle: R\sin\theta\cos\phi
		double cos_phi = vertex_normal.dot(tangent) / std::sin(theta[i]);
		if (cos_phi >=  1.0) cos_phi = 1.0;
		if (cos_phi <= -1.0) cos_phi = -1.0;
		phi[i] = std::acos(cos_phi);
		// The sine of the azimuthal angle is given as the dot product of the normal at the vertex and the
		// bitangent at the tessera center divided by the sine of the polar angle: R\sin\theta\sin\phi
		double sin_phi = vertex_normal.dot(bitangent) / std::sin(theta[i]);
		if (sin_phi <= 0.0) phi[i] = 2 * M_PI - phi[i];
		if (i != 0) {
			phi[i] -= phi[0];
			if (phi[i] < 0.0) phi[i] += 2 * M_PI;
		}
	}
	// Rewrite tangent as linear combination of original tangent and bitangent
	// then recalculate bitangent so that it's orthogonal to the tangent
	tangent = tangent * std::cos(phi[0]) + bitangent * std::sin(phi[0]);
	bitangent = normal.cross(tangent);
       
	// Populate numb and phinumb arrays
	phi[0] = 0.0;
	numb[0] = 0; numb[1] = 1;
	phinumb[0] = phi[0]; phinumb[1] = phi[1];
	for (int i = 2; i < nVertices; ++i) { // This loop is 2-based
		for (int j = 1; j < (i-1); ++j) { // This loop is 1-based 
			if (phi[i] < phinumb[j]) {
				for (int k = 0; k < (i - j); ++k) {
					numb[i - k + 1] = numb[i - k];
					phinumb[i - k + 1] = phinumb[i - k];
				}
				numb[j] = i;
				phinumb[j] = phi[i];
				goto jump; // Ugly, to be refactored!!!
			}
		}
		numb[i] = i;
		phinumb[i] = phi[i];
		jump:
		std::cout << "In definition of angles: jumped" << std::endl;
	}
	numb[nVertices] = numb[0];
	phinumb[nVertices] = 2 * M_PI;

	// Actual integration occurs here
	for (int i = 0; i < nVertices; ++i) { // Loop on edges
		// For every edge (i.e. pair of adjacent vertices) do a Gaussian product integration
		// on \phi (azimuthal integration, 64-points rule) and \theta (polar integration, 16-points rule)
		double pha = phinumb[i];
		double phb = phinumb[i+1];
		double aph = (phb - pha) / 2.0;
		double bph = (phb + pha) / 2.0;
		double tha = theta[numb[i]];
		double thb = theta[numb[i+1]];
		double thmax = 0.0;
		Eigen::Vector3d oc = (arcs.col(i) - sph.center()) / sph.radius();
		double oc_norm = oc.norm();
		double oc_norm2 = std::pow(oc_norm, 2);
		for (int j = 0; j < 32; ++j) { // Loop on Gaussian points (64-points rule)
			for (int k = 0; k <= 1; ++k) { // Obtain coordinates of Gaussian points
				double ph = (2*k - 1) * aph * gauss64Abscissa(j) + bph;
				double cos_ph = std::cos(ph);
				double sin_ph = std::sin(ph);
				if (oc_norm2 < 1.0e-07) {
					double cotg_thmax = (std::sin(ph-pha) / std::tan(thb) + std::sin(phb-ph) / std::tan(tha)) / std::sin(phb - pha);
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
				if (!(thmax < 1.0e-08)) {
					for (int l = 0; l < 8; ++l) { // Loop on Gaussian points (16-points rule)                                                         		
				        	for (int m = 0; m <= 1; ++m) { // Obtain coordinates of Gaussian points
				        		double th = (2*m - 1) * ath * gauss16Abscissa(l) + ath;
				        		double cos_th = std::cos(th);
				        		double sin_th = std::sin(th);
				        		Eigen::Vector3d point;
							point(0) = tangent(0) * sin_th * cos_ph 
				        			  + bitangent(0) * sin_th * sin_ph 
				        			  + normal(0) * (cos_th - 1.0);
				        		point(1) = tangent(1) * sin_th * cos_ph 
				        			  + bitangent(1) * sin_th * sin_ph 
				        			  + normal(1) * (cos_th - 1.0);
				        		point(2) = tangent(2) * sin_th * cos_ph 
				        			  + bitangent(2) * sin_th * sin_ph 
				        			  + normal(2) * (cos_th - 1.0);
				        		double rth = std::sqrt(2 * (1.0 - cos_th));
				        		double rtheps = gf->function(point, Eigen::Vector3d::Zero()); // Evaluate Green's function at Gaussian points
							S16 += (sph.radius() / (4 * M_PI * rtheps)) * sin_th * ath * gauss16Weight(l);
							//S16 += (sph.radius() / (4 * M_PI)) * sin_th * ath * gauss16Weight(l);
				        	}
				        }
				        S64 += S16 * aph * gauss64Weight(j);
				}
			}
		}
	}

	Sii = S64 * area;
	std::cout << "S64 = " << S64 << " area = " << area << std::endl;
	std::cout << "Sii = " << Sii << std::endl;

        return Sii;
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
	double Sii = 0.0, S64 = 0.0;
	// Extract relevant data from Element
	int nVertices = e.nVertices();
	double area = e.area();
	Eigen::Vector3d center = e.center();
	Eigen::Vector3d normal = e.normal();
	Sphere sph = e.sphere();
	Eigen::Matrix3Xd vertices = e.vertices();
	Eigen::Matrix3Xd arcs = e.arcs();

	// Calculation of the tangent and the bitangent (binormal) vectors
	Eigen::Vector3d tangent, bitangent;
	e.tangent_and_bitangent(tangent, bitangent);

	std::vector<double> theta(nVertices), phi(nVertices), phinumb(nVertices+1);
	std::vector<int> numb(nVertices+1);
	// Clean-up heap crap
	std::fill_n(theta.begin(),   nVertices,   0.0);
	std::fill_n(phi.begin(),     nVertices,   0.0);
	std::fill_n(numb.begin(),    nVertices+1, 0);
	std::fill_n(phinumb.begin(), nVertices+1, 0.0);

        // Calculate the azimuthal and polar angles for the tessera vertices:
	// we use the normal, tangent and bitangent as a local reference frame
	for (int i = 0; i < nVertices; ++i) {
		Eigen::Vector3d vertex_normal = (vertices.col(i) - sph.center()) / sph.radius();
		// The cosine of the polar angle is given as the dot product of the normal at the vertex and the
		// normal at the tessera center: R\cos\theta
		double cos_theta = vertex_normal.dot(normal);
		if (cos_theta >=  1.0) cos_theta = 1.0;
		if (cos_theta <= -1.0) cos_theta = -1.0;
		theta[i] = std::acos(cos_theta);
		// The cosine of the azimuthal angle is given as the dot product of the normal at the vertex and the
		// tangent at the tessera center divided by the sine of the polar angle: R\sin\theta\cos\phi
		double cos_phi = vertex_normal.dot(tangent) / std::sin(theta[i]);
		if (cos_phi >=  1.0) cos_phi = 1.0;
		if (cos_phi <= -1.0) cos_phi = -1.0;
		phi[i] = std::acos(cos_phi);
		// The sine of the azimuthal angle is given as the dot product of the normal at the vertex and the
		// bitangent at the tessera center divided by the sine of the polar angle: R\sin\theta\sin\phi
		double sin_phi = vertex_normal.dot(bitangent) / std::sin(theta[i]);
		if (sin_phi <= 0.0) phi[i] = 2 * M_PI - phi[i];
		if (i != 0) {
			phi[i] -= phi[0];
			if (phi[i] < 0.0) phi[i] += 2 * M_PI;
		}
	}
	// Rewrite tangent as linear combination of original tangent and bitangent
	// then recalculate bitangent so that it's orthogonal to the tangent
	tangent = tangent * std::cos(phi[0]) + bitangent * std::sin(phi[0]);
	bitangent = normal.cross(tangent);
       
	// Populate numb and phinumb arrays
	phi[0] = 0.0;
	numb[0] = 1; numb[1] = 2;
	phinumb[0] = phi[0]; phinumb[1] = phi[1];
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
		std::cout << "In definition of angles: jumped" << std::endl;
	}
	numb[nVertices] = numb[0];
	phinumb[nVertices] = 2 * M_PI;

	// Actual integration occurs here
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
					double cotg_thmax = (std::sin(ph-pha) / std::tan(thb) + std::sin(phb-ph) / std::tan(tha)) / std::sin(phb - pha);
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
				if (!(thmax < 1.0e-08)) {
					for (int l = 0; l < 8; ++l) { // Loop on Gaussian points (16-points rule)                                                         		
				        	for (int m = 0; m <= 1; ++m) {
				        		double th = (2*m - 1) * ath * gauss16Abscissa(l) + ath;
				        		double cos_th = std::cos(th);
				        		double sin_th = std::sin(th);
				        		Eigen::Vector3d point;
							point(0) = tangent(0) * sin_th * cos_ph 
				        			  + bitangent(0) * sin_th * sin_ph 
				        			  + normal(0) * (cos_th - 1.0);
				        		point(1) = tangent(1) * sin_th * cos_ph 
				        			  + bitangent(1) * sin_th * sin_ph 
				        			  + normal(1) * (cos_th - 1.0);
				        		point(2) = tangent(2) * sin_th * cos_ph 
				        			  + bitangent(2) * sin_th * sin_ph 
				        			  + normal(2) * (cos_th - 1.0);
				        		double rth = std::sqrt(2 * (1.0 - cos_th));
				        		double rtheps = gf->function(point, Eigen::Vector3d::Zero()); // Evaluate Green's function at Gaussian points
				        		S16 += (sph.radius() / (4 * M_PI * rtheps)) * sin_th * ath * gauss16Weight(l);
				        	}
				        }
				        S64 += S16 * aph * gauss64Weight(j);
				}
			}
		}
	}

	Sii = S64 * area;

        return Sii;
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
	double S64_test = 0.0;
	// Extract relevant data from Element
	int nVertices = e.nVertices();
	double area = e.area();
	Eigen::Vector3d center = e.center();
	Eigen::Vector3d normal = e.normal();
	Sphere sph = e.sphere();
	Eigen::Matrix3Xd vertices = e.vertices();
	Eigen::Matrix3Xd arcs = e.arcs();

	// Calculation of the tangent and the bitangent (binormal) vectors
	Eigen::Vector3d tangent, bitangent;
	e.tangent_and_bitangent(tangent, bitangent);

	std::vector<double> theta(nVertices), phi(nVertices), phinumb(nVertices+1);
	std::vector<int> numb(nVertices+1);
	// Clean-up heap crap
	std::fill_n(theta.begin(),   nVertices,   0.0);
	std::fill_n(phi.begin(),     nVertices,   0.0);
	std::fill_n(numb.begin(),    nVertices+1, 0);
	std::fill_n(phinumb.begin(), nVertices+1, 0.0);

        // Calculate the azimuthal and polar angles for the tessera vertices:
	// we use the normal, tangent and bitangent as a local reference frame
	for (int i = 0; i < nVertices; ++i) {
		Eigen::Vector3d vertex_normal = vertices.col(i) - sph.center();
		// The cosine of the polar angle is given as the dot product of the normal at the vertex and the
		// normal at the tessera center: R\cos\theta
		double cos_theta = vertex_normal.dot(normal) / sph.radius();
		if (cos_theta >=  1.0) cos_theta = 1.0;
		if (cos_theta <= -1.0) cos_theta = -1.0;
		theta[i] = std::acos(cos_theta);
		// The cosine of the azimuthal angle is given as the dot product of the normal at the vertex and the
		// tangent at the tessera center divided by the sine of the polar angle: R\sin\theta\cos\phi
		double cos_phi = vertex_normal.dot(tangent) / (sph.radius() * std::sin(theta[i]));
		if (cos_phi >=  1.0) cos_phi = 1.0;
		if (cos_phi <= -1.0) cos_phi = -1.0;
		phi[i] = std::acos(cos_phi);
		// The sine of the azimuthal angle is given as the dot product of the normal at the vertex and the
		// bitangent at the tessera center divided by the sine of the polar angle: R\sin\theta\sin\phi
		double sin_phi = vertex_normal.dot(bitangent) / (sph.radius() * std::sin(theta[i]));
		if (sin_phi <= 0.0) phi[i] = 2 * M_PI - phi[i];
	}
	for (int i = 1; i < nVertices; ++i) {
		phi[i] = phi[i] - phi[0];
		if (phi[i] < 0.0) phi[i] = 2 * M_PI + phi[i];
	}
	// Rewrite tangent as linear combination of original tangent and bitangent
	// then recalculate bitangent so that it's orthogonal to the tangent
	tangent = tangent * std::cos(phi[0]) + bitangent * std::sin(phi[0]);
	bitangent = normal.cross(tangent);
       
	// Populate numb and phinumb arrays
	phi[0] = 0.0;
	numb[0] = 0; numb[1] = 1;
	phinumb[0] = phi[0]; phinumb[1] = phi[1];
	for (int i = 2; i < nVertices; ++i) { // This loop is 2-based
		for (int j = 1; j < i; ++j) { // This loop is 1-based 
			if (phi[i] < phinumb[j]) {
				for (int k = 0; k < (i - j); ++k) {
					numb[i - k] = numb[i - k -1];
					phinumb[i - k] = phinumb[i - k -1];
				}
				numb[j] = i;
				phinumb[j] = phi[i];
				goto jump; // Ugly, to be refactored!!!
			}
		}
		numb[i] = i;
		phinumb[i] = phi[i];
		jump:
		std::cout << "In definition of angles: jumped" << std::endl;
	}
	numb[nVertices] = numb[0];
	phinumb[nVertices] = 2 * M_PI;
	
	// Actual integration occurs here
	for (int i = 0; i < nVertices; ++i) { // Loop on edges
		double pha = phinumb[i];
		double phb = phinumb[i+1];
		double aph = (phb - pha) / 2.0;
		double bph = (phb + pha) / 2.0;
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
					double cotg_thmax = (std::sin(ph-pha) / std::tan(thb) + std::sin(phb-ph) / std::tan(tha)) / std::sin(phb - pha);
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
				double S16_test = 0.0;
				if (!(thmax < 1.0e-08)) {
					for (int l = 0; l < 8; ++l) { // Loop on Gaussian points (16-points rule)                                                         		
				        	for (int m = 0; m <= 1; ++m) {
				        		double th = (2*m - 1) * ath * gauss16Abscissa(l) + ath;
				        		double cos_th = std::cos(th);
				        		double sin_th = std::sin(th);
				        		Eigen::Vector3d point;
							point(0) = tangent(0) * sin_th * cos_ph 
				        			  + bitangent(0) * sin_th * sin_ph 
				        			  + normal(0) * (cos_th - 1.0);
				        		point(1) = tangent(1) * sin_th * cos_ph 
				        			  + bitangent(1) * sin_th * sin_ph 
				        			  + normal(1) * (cos_th - 1.0);
				        		point(2) = tangent(2) * sin_th * cos_ph 
				        			  + bitangent(2) * sin_th * sin_ph 
				        			  + normal(2) * (cos_th - 1.0);
				        		double rth = std::sqrt(2 * (1.0 - cos_th));
				        		double rtheps = gf->function(point, Eigen::Vector3d::Zero()); // Evaluate Green's function at Gaussian points
				        		S16 += (sph.radius() / rtheps) * sin_th * ath * gauss16Weight(l);
							S16_test += sph.radius()*sin_th*ath*gauss16Weight(l);
				        	}
				        }
				        S64 += S16 * aph * gauss64Weight(j);
					S64_test += S16_test * aph * gauss64Weight(j);
				}
			}
		}
	}

	Sii = S64;
        std::cout << "e.area() = " << area << std::endl;
	std::cout << "S64_test = " << S64_test << std::endl;
	std::cout << "diff = " << area-S64_test << std::endl;

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
