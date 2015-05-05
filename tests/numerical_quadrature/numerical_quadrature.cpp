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

#define BOOST_TEST_MODULE NumericalQuadrature

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <cmath>
#include <iomanip>
#include <limits>

#include "Config.hpp"

#include <Eigen/Dense>

#include <boost/bind.hpp>
#include <boost/bind/placeholders.hpp>
#include <boost/function.hpp>

#include "cnpyPimpl.hpp"
#include "Element.hpp"
#include "GePolCavity.hpp"
#include "MathUtils.hpp"
#include "PhysicalConstants.hpp"
#include "TestingMolecules.hpp"

typedef boost::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)>
singleLayerIntegrand;
typedef boost::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &,
		        const Eigen::Vector3d &)> doubleLayerIntegrand;

/*! \class NumericalQuadrature
 *  \test \b NumericalQuadrature_sphere tests numerical quadrature on a sphere for integrand = 1.0
 */
BOOST_AUTO_TEST_CASE(sphere)
{
    struct f {
	double operator()(const Eigen::Vector3d & /* s */, const Eigen::Vector3d & /* p */) { return 1.0; }
    };

    double radius = 1.55;
    double area = 0.4;
    Molecule point = dummy<0>(radius);
    GePolCavity cavity(point, area, 0.0, 100.0);
    Eigen::VectorXd results = Eigen::VectorXd::Zero(cavity.size());

    singleLayerIntegrand F = f();

    for (int i = 0; i < cavity.size(); ++i) {
	results(i) = integrator<32, 16>(F, cavity.elements(i));
	double diff = results(i) - cavity.elementArea(i);
	if (std::abs(diff) > 1.0e-12) {
	    BOOST_TEST_MESSAGE("Test versus area for single sphere");
            BOOST_TEST_MESSAGE("Tessera n. " << i+1);
	    BOOST_TEST_MESSAGE("diff = " << results(i) - cavity.elementArea(i));
	}
    }

    for (int i = 0; i < cavity.size(); ++i) {
    	BOOST_REQUIRE_CLOSE(results(i), cavity.elementArea(i), 1.0e-12);
    }
}

/*! \class NumericalQuadrature
 *  \test \b NumericalQuadrature_sphere_1r tests numerical quadrature on a sphere for integrand 1.0/r
 */
BOOST_AUTO_TEST_CASE(sphere_1r)
{
    double radius = 1.55;
    double area = 0.4;
    Molecule point = dummy<0>(radius);
    GePolCavity cavity(point, area, 0.0, 100.0);

    struct f {
	double operator()(const Eigen::Vector3d & /* s */, const Eigen::Vector3d & /* p */) { double radius = 1.55; return (1.0/radius); }
    };

    Eigen::VectorXd results = Eigen::VectorXd::Zero(cavity.size());

    singleLayerIntegrand F = f();

    for (int i = 0; i < cavity.size(); ++i) {
	results(i) = integrator<32, 16>(F, cavity.elements(i));
	double diff = results(i) - (cavity.elementArea(i)/radius);
	if (std::abs(diff) > 1.0e-12) {
	    BOOST_TEST_MESSAGE("Test versus area divided by radius for single sphere");
            BOOST_TEST_MESSAGE("Tessera n. " << i+1);
	    BOOST_TEST_MESSAGE("diff = " << results(i) - (cavity.elementArea(i)/radius) );
	}
    }

    for (int i = 0; i < cavity.size(); ++i) {
    	BOOST_REQUIRE_CLOSE(results(i), (cavity.elementArea(i)/radius), 1.0e-12);
    }
}

/*! \class NumericalQuadrature
 *  \test \b NumericalQuadrature_molecule tests numerical quadrature function on H2 cavity
 */
BOOST_AUTO_TEST_CASE(molecule)
{
    struct f {
	double operator()(const Eigen::Vector3d & /* s */, const Eigen::Vector3d & /* p */) { return 1.0; }
    };

    double area = 0.2;
    double probeRadius = 1.385;
    double minRadius = 0.2;
    Molecule molec = H2();
    GePolCavity cavity(molec, area, probeRadius, minRadius);
    Eigen::VectorXd results = Eigen::VectorXd::Zero(cavity.size());

    singleLayerIntegrand F = f();

    for (int i = 0; i < cavity.size(); ++i) {
	results(i) = integrator<64, 16>(F, cavity.elements(i));
	double diff = results(i) - cavity.elementArea(i);
	if (std::abs(diff) > 1.0e-12) {
	    BOOST_TEST_MESSAGE("Test versus area for H2 molecule");
            BOOST_TEST_MESSAGE("Tessera n. " << i+1);
	    BOOST_TEST_MESSAGE("diff = " << results(i) - cavity.elementArea(i));
	}
    }
    /*
    // In case you need to update the reference files...
    unsigned int dim = static_cast<unsigned int>(cavity.size());
    const unsigned int shape[] = {dim};
    cnpy::npy_save("molecule.npy", results.data(), shape, 1, "w", false);
    */
    cnpy::NpyArray raw_ref = cnpy::npy_load("molecule.npy");
    int dim = raw_ref.shape[0];
    Eigen::VectorXd reference = Eigen::VectorXd::Zero(dim);
    reference = getFromRawBuffer<double>(dim, 1, raw_ref.data);

    for (int i = 0; i < cavity.size(); ++i) {
    	BOOST_REQUIRE_CLOSE(results(i), reference(i), 1.0e-12);
    }
}

/*! \class NumericalQuadrature
 *  \test \b NumericalQuadrature_molecule_1r tests numerical quadrature function on H2 cavity
 */
BOOST_AUTO_TEST_CASE(molecule_1r)
{
    struct f {
	double operator()(const Eigen::Vector3d & /* s */, const Eigen::Vector3d & /* p */) { double radius = 1.20; return (1.0/radius); }
    };

    double area = 0.2;
    double probeRadius = 1.385;
    double minRadius = 0.2;
    Molecule molec = H2();
    Sphere sph1 = molec.spheres(0);
    GePolCavity cavity(molec, area, probeRadius, minRadius);
    Eigen::VectorXd results = Eigen::VectorXd::Zero(cavity.size());

    singleLayerIntegrand F = f();

    for (int i = 0; i < cavity.size(); ++i) {
	results(i) = integrator<64, 16>(F, cavity.elements(i));
	double diff = results(i) - (cavity.elementArea(i)/sph1.radius());
	if (std::abs(diff) > 1.0e-12) {
	    BOOST_TEST_MESSAGE("Test versus area divided by radius for H2 molecule");
            BOOST_TEST_MESSAGE("Tessera n. " << i+1);
	    BOOST_TEST_MESSAGE("diff = " << results(i) - (cavity.elementArea(i)/sph1.radius()) );
	}
    }
    /*
    // In case you need to update the reference files...
    unsigned int dim = static_cast<unsigned int>(cavity.size());
    const unsigned int shape[] = {dim};
    cnpy::npy_save("molecule_1r.npy", results.data(), shape, 1, "w", false);
    */
    cnpy::NpyArray raw_ref = cnpy::npy_load("molecule_1r.npy");
    int dim = raw_ref.shape[0];
    Eigen::VectorXd reference = Eigen::VectorXd::Zero(dim);
    reference = getFromRawBuffer<double>(dim, 1, raw_ref.data);

    for (int i = 0; i < cavity.size(); ++i) {
    	BOOST_REQUIRE_CLOSE(results(i), reference(i), 1.0e-12);
    }
}
