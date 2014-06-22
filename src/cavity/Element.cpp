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

#include "Element.hpp"

#include <cmath>
#include <stdexcept>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include <boost/math/special_functions/sign.hpp>

#include "Sphere.hpp"

void Element::tangent_and_bitangent(Eigen::Vector3d & t_, Eigen::Vector3d & b_) const
{
        double rmin = 0.99;
        if (std::fabs(normal_(0)) <= rmin) {
        	rmin = std::fabs(normal_(0));
        	t_(0) = 0.0;
        	t_(1) = - normal_(2) / std::sqrt(1.0 - std::pow(normal_(0), 2));
        	t_(2) = normal_(1) / std::sqrt(1.0 - std::pow(normal_(0), 2));
        }
        if (std::fabs(normal_(1)) <= rmin) {
        	rmin = std::fabs(normal_(1));
        	t_(0) = normal_(2) / std::sqrt(1.0 - std::pow(normal_(1), 2));
        	t_(1) = 0.0;
        	t_(2) = -normal_(0) / std::sqrt(1.0 - std::pow(normal_(1), 2));
        }
        if (std::fabs(normal_(2)) <= rmin) {
        	rmin = std::fabs(normal_(2));
        	t_(0) = normal_(1) / std::sqrt(1.0 - std::pow(normal_(2), 2));
        	t_(1) = -normal_(0) / std::sqrt(1.0 - std::pow(normal_(2), 2));
        	t_(2) =0.0; 
        }
        b_ = normal_.cross(t_);
	// Check that the calculated Frenet-Serret frame is left-handed (levogiro)
	// by checking that the determinant of the matrix whose columns are the normal,
	// tangent and bitangent vectors has determinant 1 (the system is orthonormal!)
	Eigen::Matrix3d M;
        M.col(0) = normal_;
	M.col(1) = t_;
	M.col(2) = b_;
	if (boost::math::sign(M.determinant()) != 1) {
		throw std::runtime_error("Frenet-Serret local frame is not left-handed!");
	}
}
