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

#ifndef MATHUTILS_HPP
#define MATHUTILS_HPP

#include <algorithm>
#include <cmath>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include <boost/bind.hpp>
#include <boost/bind/placeholders.hpp>
#include <boost/function.hpp>

#include "Element.hpp"
#include "Sphere.hpp"
#include "Symmetry.hpp"

/*! \fn inline bool isZero(double value, double threshold)
 *  \param[in] value     the value to be checked
 *  \param[in] threshold the threshold
 *
 *  Returns true if value is less or equal to threshold
 */
inline bool isZero(double value, double threshold)
{
    return (std::abs(value) <= threshold);
}

/*! \fn inline bool numericalZero(double value)
 *  \param[in] value the value to be checked
 *
 *  Returns true if value is less than 1.0e-14
 */
inline bool numericalZero(double value)
{
    return (isZero(value, 1.0e-14));
}

/*! \fn inline void symmetryBlocking(Eigen::MatrixXd & matrix, int cavitySize, int ntsirr, int nr_irrep)
 *  \param[out] matrix the matrix to be block-diagonalized
 *  \param[in]  cavitySize the size of the cavity (size of the matrix)
 *  \param[in]  ntsirr     the size of the irreducible portion of the cavity (size of the blocks)
 *  \param[in]  nr_irrep   the number of irreducible representations (number of blocks)
 */
inline void symmetryBlocking(Eigen::MatrixXd & matrix, int cavitySize, int ntsirr,
                             int nr_irrep)
{
    // This function implements the simmetry-blocking of the PCM
    // matrix due to point group symmetry as reported in:
    // L. Frediani, R. Cammi, C. S. Pomelli, J. Tomasi and K. Ruud, J. Comput.Chem. 25, 375 (2003)
    // u is the character table for the group (t in the paper)
    Eigen::MatrixXd u = Eigen::MatrixXd::Zero(nr_irrep, nr_irrep);
    for (int i = 0; i < nr_irrep; ++i) {
        for (int j = 0; j < nr_irrep; ++j) {
            u(i, j) = Symmetry::parity(i&j);
        }
    }
    // Naming of indices:
    //     a, b, c, d   run over the total size of the cavity (N)
    //     i, j, k, l   run over the number of irreps (n)
    //     p, q, r, s   run over the irreducible size of the cavity (N/n)
    // Instead of forming U (T in the paper) and then perform the dense
    // multiplication, we multiply block-by-block using just the u matrix.
    //      matrix = U * matrix * Ut; U * Ut = Ut * U = id
    // First half-transformation, i.e. first_half = matrix * Ut
    Eigen::MatrixXd first_half = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
    for (int i = 0; i < nr_irrep; ++i) {
        int ioff = i * ntsirr;
        for (int k = 0; k < nr_irrep; ++k) {
            int koff = k * ntsirr;
            for (int j = 0; j < nr_irrep; ++j) {
                int joff = j * ntsirr;
                double ujk = u(j, k) / nr_irrep;
                for (int p = 0; p < ntsirr; ++p) {
                    int a = ioff + p;
                    for (int q = 0; q < ntsirr; ++q) {
                        int b = joff + q;
                        int c = koff + q;
                        first_half(a, c) += matrix(a, b) * ujk;
                    }
                }
            }
        }
    }
    // Second half-transformation, i.e. matrix = U * first_half
    matrix.setZero(cavitySize, cavitySize);
    for (int i = 0; i < nr_irrep; ++i) {
        int ioff = i * ntsirr;
        for (int k = 0; k < nr_irrep; ++k) {
            int koff = k * ntsirr;
            for (int j = 0; j < nr_irrep; ++j) {
                int joff = j * ntsirr;
                double uij = u(i, j);
                for (int p = 0; p < ntsirr; ++p) {
                    int a = ioff + p;
                    int b = joff + p;
                    for (int q = 0; q < ntsirr; ++q) {
                        int c = koff + q;
                        matrix(a, c) += uij * first_half(b, c);
                    }
                }
            }
        }
    }
    // Traverse the matrix and discard numerical zeros
    for (int a = 0; a < cavitySize; ++a) {
        for (int b = 0; b < cavitySize; ++b) {
            if (numericalZero(matrix(a, b))) {
                matrix(a, b) = 0.0;
            }
        }
    }
}

/*! \fn inline void symmetryPacking(std::vector<Eigen::MatrixXd> & blockedMatrix, const Eigen::MatrixXd & fullMatrix, int nrBlocks, int dimBlock)
 *  \param[out] blockedMatrix the result of packing fullMatrix
 *  \param[in]  fullMatrix the matrix to be packed
 *  \param[in]  dimBlock the dimension of the square blocks
 *  \param[in]  nrBlocks the number of square blocks
 */
inline void symmetryPacking(std::vector<Eigen::MatrixXd> & blockedMatrix,
                            const Eigen::MatrixXd & fullMatrix, int dimBlock, int nrBlocks)
{
    // This function packs the square block diagonal fullMatrix with nrBlocks of dimension dimBlock
    // into a std::vector<Eigen::MatrixXd> containing nrBlocks square matrices of dimenion dimBlock.
    int j = 0;
    for (int i = 0; i < nrBlocks; ++i) {
        blockedMatrix.push_back(fullMatrix.block(j, j, dimBlock, dimBlock));
        j += dimBlock;
    }
}

/*! \fn inline void hermitivitize(Eigen::MatrixBase<Derived> & matrix_)
 *  \param[out] matrix_ the matrix to be hermitivitized
 *  \tparam     Derived the numeric type of matrix_ elements
 *
 *  Given matrix_ returns 0.5 * (matrix_ + matrix_^dagger)
 */
template <typename Derived>
inline void hermitivitize(Eigen::MatrixBase<Derived> & matrix_)
{
    // We need to use adjoint().eval() to avoid aliasing issues, see:
    // http://eigen.tuxfamily.org/dox/group__TopicAliasing.html
    // The adjoint is evaluated explicitly into an intermediate.
    matrix_ = 0.5 * (matrix_ + matrix_.adjoint().eval());
}

/*! \fn inline void eulerRotation(Eigen::Matrix3d & R_, const Eigen::Vector3d & eulerAngles_)
 *  \brief Build rotation matrix between two reference frames given the Euler angles.
 *  \param[out] R_ the rotation matrix
 *  \param[in]  eulerAngles_ the Euler angles, in degrees, describing the rotation
 *
 *  We assume the convention \f$ R = Z_3 X_2 Z_1 \f$ for the ordering of the extrinsic 
 *  elemental rotations (see http://en.wikipedia.org/wiki/Euler_angles)
 *  The Euler angles are given in the order \f$ \phi, \theta, \psi \f$. 
 *  If we write \f$ c_i, s_i \,\, i = 1, 3 \f$ for their cosines and sines the rotation
 *  matrix will be:
 *  \f[
 *  	R = \begin{pmatrix}
 *  	      c_1c_3 - s_1c_2s_3 & -s_1c_3 - c_1c_2s_3 &  s_2s_3 \\
 *  	      c_1s_3 + s_1c_2c_3 & -s_1s_3 + c_1c_2c_3 & -s_2c_3 \\
 *  	      s_1s_2             & c_1s_2              &  c_2
 *  	    \end{pmatrix}
 *  \f]
 *  Eigen's geometry module is used to calculate the rotation matrix
 */
inline void eulerRotation(Eigen::Matrix3d & R_, const Eigen::Vector3d & eulerAngles_)
{
	double to_radians = M_PI / 180.0;
	double phi   = eulerAngles_(0) * to_radians;
	double theta = eulerAngles_(1) * to_radians;
	double psi   = eulerAngles_(2) * to_radians;
	R_ = Eigen::AngleAxis<double>(psi,   Eigen::Vector3d::UnitZ())
           * Eigen::AngleAxis<double>(theta, Eigen::Vector3d::UnitX())
           * Eigen::AngleAxis<double>(phi,   Eigen::Vector3d::UnitZ());
}

/*! Abscissae for 16-point Gaussian quadrature rule */
inline double gauss16Abscissa(int i)
{
   static std::vector<double> x16(8);

   x16[0] = 0.9894009349916499325961542;
   x16[1] = 0.9445750230732325760779884;
   x16[2] = 0.8656312023878317438804679;
   x16[3] = 0.7554044083550030338951012;
   x16[4] = 0.6178762444026437484466718;
   x16[5] = 0.4580167776572273863424194;
   x16[6] = 0.2816035507792589132304605;
   x16[7] = 0.0950125098376374401853193;

   return x16[i];
}

/*! Weights for 16-point Gaussian quadrature rule */
inline double gauss16Weight(int i)
{
   static std::vector<double> w16(8);

   w16[0] = 0.0271524594117540948517806;
   w16[1] = 0.0622535239386478928628438;
   w16[2] = 0.0951585116824927848099251;
   w16[3] = 0.1246289712555338720524763;
   w16[4] = 0.1495959888165767320815017;
   w16[5] = 0.1691565193950025381893121;
   w16[6] = 0.1826034150449235888667637;
   w16[7] = 0.1894506104550684962853967;

   return w16[i];
}

/*! Abscissae for 64-point Gaussian quadrature rule */
inline double gauss64Abscissa(int i)
{
   static std::vector<double> x64(32);

   x64[0]  = 0.9993050417357721394569056; 
   x64[1]  = 0.9963401167719552793469245;
   x64[2]  = 0.9910133714767443207393824;
   x64[3]  = 0.9833362538846259569312993;
   x64[4]  = 0.9733268277899109637418535;
   x64[5]  = 0.9610087996520537189186141;
   x64[6]  = 0.9464113748584028160624815;
   x64[7]  = 0.9295691721319395758214902;
   x64[8]  = 0.9105221370785028057563807;
   x64[9]  = 0.8893154459951141058534040;
   x64[10] = 0.8659993981540928197607834;
   x64[11] = 0.8406292962525803627516915;
   x64[12] = 0.8132653151227975597419233;
   x64[13] = 0.7839723589433414076102205;
   x64[14] = 0.7528199072605318966118638;
   x64[15] = 0.7198818501716108268489402;
   x64[16] = 0.6852363130542332425635584;
   x64[17] = 0.6489654712546573398577612;
   x64[18] = 0.6111553551723932502488530;
   x64[19] = 0.5718956462026340342838781;
   x64[20] = 0.5312794640198945456580139;
   x64[21] = 0.4894031457070529574785263;
   x64[22] = 0.4463660172534640879849477;
   x64[23] = 0.4022701579639916036957668;
   x64[24] = 0.3572201583376681159504426;
   x64[25] = 0.3113228719902109561575127;
   x64[26] = 0.2646871622087674163739642;
   x64[27] = 0.2174236437400070841496487;
   x64[28] = 0.1696444204239928180373136;
   x64[29] = 0.1214628192961205544703765;
   x64[30] = 0.0729931217877990394495429;
   x64[31] = 0.0243502926634244325089558;

   return x64[i];
}

/*! Weights for 64-point Gaussian quadrature rule */
inline double gauss64Weight(int i)
{
   static std::vector<double> w64(32);

   w64[0]  = 0.0017832807216964329472961; 
   w64[1]  = 0.0041470332605624676352875;
   w64[2]  = 0.0065044579689783628561174;
   w64[3]  = 0.0088467598263639477230309;
   w64[4]  = 0.0111681394601311288185905;
   w64[5]  = 0.0134630478967186425980608;
   w64[6]  = 0.0157260304760247193219660;
   w64[7]  = 0.0179517157756973430850453;
   w64[8]  = 0.0201348231535302093723403;
   w64[9]  = 0.0222701738083832541592983;
   w64[10] = 0.0243527025687108733381776;
   w64[11] = 0.0263774697150546586716918;
   w64[12] = 0.0283396726142594832275113;
   w64[13] = 0.0302346570724024788679741;
   w64[14] = 0.0320579283548515535854675;
   w64[15] = 0.0338051618371416093915655;
   w64[16] = 0.0354722132568823838106931;
   w64[17] = 0.0370551285402400460404151;
   w64[18] = 0.0385501531786156291289625;
   w64[19] = 0.0399537411327203413866569;
   w64[20] = 0.0412625632426235286101563;
   w64[21] = 0.0424735151236535890073398;
   w64[22] = 0.0435837245293234533768279;
   w64[23] = 0.0445905581637565630601347;
   w64[24] = 0.0454916279274181444797710;
   w64[25] = 0.0462847965813144172959532;
   w64[26] = 0.0469681828162100173253263;
   w64[27] = 0.0475401657148303086622822;
   w64[28] = 0.0479993885964583077281262;
   w64[29] = 0.0483447622348029571697695;
   w64[30] = 0.0485754674415034269347991;
   w64[31] = 0.0486909570091397203833654;

   return w64[i];
}

typedef boost::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)> integrand;

inline double integrator(const integrand & F_, const Element & e)
{
	double S64 = 0.0;
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
				        		double rtheps = F_(point, Eigen::Vector3d::Zero()); // Evaluate Green's function at Gaussian points
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

	std::cout << "In integrator!!" << std::endl;
        std::cout << "e.area() = " << area << std::endl;
	std::cout << "S64_test = " << S64_test << std::endl;
	std::cout << "diff = " << area-S64_test << std::endl;

	return S64;
}

#endif // MATHUTILS_HPP
