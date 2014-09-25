#define BOOST_TEST_MODULE NumericalIntegrator 

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "cnpyPimpl.hpp"
#include "DerivativeTypes.hpp"
#include "AnisotropicLiquid.hpp"
#include "Element.hpp"
#include "GePolCavity.hpp"
#include "IonicLiquid.hpp"
#include "NumericalIntegrator.hpp"
#include "PhysicalConstants.hpp"
#include "UniformDielectric.hpp"
#include "Vacuum.hpp"

struct NumericalIntegratorTest {
    double radius;
    Eigen::Vector3d epsilon;
    Eigen::Vector3d euler;
    double eps;
    double kappa;
    Eigen::Vector3d source, probe, sourceNormal, probeNormal;
    GePolCavity cavity;
    NumericalIntegrator * diag;
    NumericalIntegratorTest() { SetUp(); }
    void SetUp() {
/*    	epsilon << 2.0, 80.0, 15.0;
    	euler << 6.0, 40.0, 15.0;*/
	epsilon << 80.0, 80.0, 80.0;
	euler << 0.0, 0.0, 0.0;
	eps = 80.0;
	kappa = 1.5;
        source = Eigen::Vector3d::Random();
        sourceNormal = source + Eigen::Vector3d::Random();
        sourceNormal.normalize();
        probe = Eigen::Vector3d::Random();
        probeNormal = probe + Eigen::Vector3d::Random();
        probeNormal.normalize();
        
	Eigen::Vector3d origin(0.0, 0.0, 0.0);
        std::vector<Sphere> spheres;
	radius = 1.44 / convertBohrToAngstrom;
        Sphere sph1(origin,  radius);
        spheres.push_back(sph1);
        double area = 10.0;
        // C1
        Symmetry pGroup = buildGroup(0, 0, 0, 0);
        cavity = GePolCavity(spheres, area, 0.0, 100.0, pGroup);
	
	diag = new NumericalIntegrator();
    }
};

/*! \class NumericalIntegrator 
 *  \test \b NumericalIntegratorTest_vacuum tests the numerical evaluation of the vacuum diagonal elements of S and D
 */
BOOST_FIXTURE_TEST_CASE(vacuum, NumericalIntegratorTest)
{
    Vacuum<AD_directional> gf(diag);

    BOOST_TEST_MESSAGE("Vacuum");
    Eigen::VectorXd S_results = Eigen::VectorXd::Zero(cavity.size());
    for (int i = 0; i < cavity.size(); ++i) {
	S_results(i) = gf.diagonalS(cavity.elements(i)); 
    }
    // In case you need to update the reference files...
    /*
    unsigned int dim = static_cast<unsigned int>(cavity.size());
    const unsigned int shape[] = {dim};
    cnpy::npy_save("vacuum_S_numerical.npy", S_results.data(), shape, 1, "w", false);
    cnpy::NpyArray raw_S_ref = cnpy::npy_load("vacuum_S_numerical.npy");
    int dim = raw_S_ref.shape[0];
    Eigen::VectorXd S_reference = Eigen::VectorXd::Zero(dim);
    S_reference = cnpy::getFromRawBuffer<double>(dim, 1, raw_S_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
    	BOOST_REQUIRE_CLOSE(S_results(i), S_reference(i), 1.0e-12);
    }
    */
    BOOST_TEST_MESSAGE("S operator diagonal by numerical quadrature");
    for (int i = 0; i < cavity.size(); ++i) {
	    BOOST_TEST_MESSAGE("S_{" << i+1 << ", " << i+1 << "} = " 
			    << std::setprecision(std::numeric_limits<long double>::digits10) << S_results(i));
    }
    
    Eigen::VectorXd D_results = Eigen::VectorXd::Zero(cavity.size());
    for (int i = 0; i < cavity.size(); ++i) {
	D_results(i) = gf.diagonalD(cavity.elements(i)); 
    }
    // In case you need to update the reference files...
    /*
    cnpy::npy_save("vacuum_D_numerical.npy", D_results.data(), shape, 1, "w", false);
    cnpy::NpyArray raw_D_ref = cnpy::npy_load("vacuum_D_numerical.npy");
    dim = raw_D_ref.shape[0];
    Eigen::VectorXd D_reference = Eigen::VectorXd::Zero(dim);
    D_reference = cnpy::getFromRawBuffer<double>(dim, 1, raw_D_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
    	BOOST_REQUIRE_CLOSE(D_results(i), D_reference(i), 1.0e-12);
    }
    */
    BOOST_TEST_MESSAGE("D operator diagonal by numerical quadrature");
    for (int i = 0; i < cavity.size(); ++i) {
	    BOOST_TEST_MESSAGE("D_{" << i+1 << ", " << i+1 << "} = " 
			    << std::setprecision(std::numeric_limits<long double>::digits10) << D_results(i));
    }
}

/*! \class NumericalIntegrator 
 *  \test \b NumericalIntegratorTest_uniformdielectric tests the numerical evaluation of the uniform dielectric diagonal elements of S and D
 */
BOOST_FIXTURE_TEST_CASE(uniformdielectric, NumericalIntegratorTest)
{
    UniformDielectric<AD_directional> gf(eps, diag);
    
    BOOST_TEST_MESSAGE("UniformDielectric");
    Eigen::VectorXd S_results = Eigen::VectorXd::Zero(cavity.size());
    for (int i = 0; i < cavity.size(); ++i) {
	S_results(i) = gf.diagonalS(cavity.elements(i)); 
    }
    // In case you need to update the reference files...
    /*
    unsigned int dim = static_cast<unsigned int>(cavity.size());
    const unsigned int shape[] = {dim};
    cnpy::npy_save("uniformdielectric_S_numerical.npy", S_results.data(), shape, 1, "w", false);
    cnpy::NpyArray raw_S_ref = cnpy::npy_load("uniformdielectric_S_numerical.npy");
    int dim = raw_S_ref.shape[0];
    Eigen::VectorXd S_reference = Eigen::VectorXd::Zero(dim);
    S_reference = cnpy::getFromRawBuffer<double>(dim, 1, raw_S_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
    	BOOST_REQUIRE_CLOSE(S_results(i), S_reference(i), 1.0e-12);
    }
    */
    BOOST_TEST_MESSAGE("S operator diagonal by numerical quadrature");
    for (int i = 0; i < cavity.size(); ++i) {
	    BOOST_TEST_MESSAGE("S_{" << i+1 << ", " << i+1 << "} = " 
			    << std::setprecision(std::numeric_limits<long double>::digits10) << S_results(i));
    }
    
    Eigen::VectorXd D_results = Eigen::VectorXd::Zero(cavity.size());
    for (int i = 0; i < cavity.size(); ++i) {
	D_results(i) = gf.diagonalD(cavity.elements(i)); 
    }
    // In case you need to update the reference files...
    /*
    cnpy::npy_save("uniformdielectric_D_numerical.npy", D_results.data(), shape, 1, "w", false);
    cnpy::NpyArray raw_D_ref = cnpy::npy_load("uniformdielectric_D_numerical.npy");
    dim = raw_D_ref.shape[0];
    Eigen::VectorXd D_reference = Eigen::VectorXd::Zero(dim);
    D_reference = cnpy::getFromRawBuffer<double>(dim, 1, raw_D_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
    	BOOST_REQUIRE_CLOSE(D_results(i), D_reference(i), 1.0e-12);
    }
    */
    BOOST_TEST_MESSAGE("D operator diagonal by numerical quadrature");
    for (int i = 0; i < cavity.size(); ++i) {
	    BOOST_TEST_MESSAGE("D_{" << i+1 << ", " << i+1 << "} = "
			    << std::setprecision(std::numeric_limits<long double>::digits10) << D_results(i));
    }
}

/*! \class NumericalIntegrator 
 *  \test \b NumericalIntegratorTest_ionic tests the numerical evaluation of the ionic liquid diagonal elements of S and D
 */
BOOST_FIXTURE_TEST_CASE(ionic, NumericalIntegratorTest)
{
    IonicLiquid<AD_directional> gf(eps, kappa, diag);

    BOOST_TEST_MESSAGE("IonicLiquid");
    Eigen::VectorXd S_results = Eigen::VectorXd::Zero(cavity.size());
    for (int i = 0; i < cavity.size(); ++i) {
	S_results(i) = gf.diagonalS(cavity.elements(i)); 
    }
    // In case you need to update the reference files...
    /*
    unsigned int dim = static_cast<unsigned int>(cavity.size());
    const unsigned int shape[] = {dim};
    cnpy::npy_save("ionic_S_numerical.npy", S_results.data(), shape, 1, "w", false);
    cnpy::NpyArray raw_S_ref = cnpy::npy_load("ionic_S_numerical.npy");
    int dim = raw_S_ref.shape[0];
    Eigen::VectorXd S_reference = Eigen::VectorXd::Zero(dim);
    S_reference = cnpy::getFromRawBuffer<double>(dim, 1, raw_S_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
    	BOOST_REQUIRE_CLOSE(S_results(i), S_reference(i), 1.0e-12);
    }
    */
    BOOST_TEST_MESSAGE("S operator diagonal by numerical quadrature");
    for (int i = 0; i < cavity.size(); ++i) {
	    BOOST_TEST_MESSAGE("S_{" << i+1 << ", " << i+1 << "} = " 
			    << std::setprecision(std::numeric_limits<long double>::digits10) << S_results(i));
    }
    
    Eigen::VectorXd D_results = Eigen::VectorXd::Zero(cavity.size());
    for (int i = 0; i < cavity.size(); ++i) {
	D_results(i) = gf.diagonalD(cavity.elements(i)); 
    }
    // In case you need to update the reference files...
    /*
    cnpy::npy_save("ionic_D_numerical.npy", D_results.data(), shape, 1, "w", false);
    cnpy::NpyArray raw_D_ref = cnpy::npy_load("ionic_D_numerical.npy");
    dim = raw_D_ref.shape[0];
    Eigen::VectorXd D_reference = Eigen::VectorXd::Zero(dim);
    D_reference = cnpy::getFromRawBuffer<double>(dim, 1, raw_D_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
    	BOOST_REQUIRE_CLOSE(D_results(i), D_reference(i), 1.0e-12);
    }
    */
    BOOST_TEST_MESSAGE("D operator diagonal by numerical quadrature");
    for (int i = 0; i < cavity.size(); ++i) {
	    BOOST_TEST_MESSAGE("D_{" << i+1 << ", " << i+1 << "} = "
			    << std::setprecision(std::numeric_limits<long double>::digits10) << D_results(i));
    }
}

/*! \class NumericalIntegrator 
 *  \test \b NumericalIntegratorTest_anisotropic tests the numerical evaluation of the anisotropic liquid diagonal elements of S and D
 */
BOOST_FIXTURE_TEST_CASE(anisotropic, NumericalIntegratorTest)
{
    AnisotropicLiquid<AD_directional> gf(epsilon, euler, diag);

    BOOST_TEST_MESSAGE("AnisotropicLiquid");
    Eigen::VectorXd S_results = Eigen::VectorXd::Zero(cavity.size());
    for (int i = 0; i < cavity.size(); ++i) {
	S_results(i) = gf.diagonalS(cavity.elements(i)); 
    }
    // In case you need to update the reference files...
    /*
    unsigned int dim = static_cast<unsigned int>(cavity.size());
    const unsigned int shape[] = {dim};
    cnpy::npy_save("anisotropic_S_numerical.npy", S_results.data(), shape, 1, "w", false);
    cnpy::NpyArray raw_S_ref = cnpy::npy_load("anisotropic_S_numerical.npy");
    int dim = raw_S_ref.shape[0];
    Eigen::VectorXd S_reference = Eigen::VectorXd::Zero(dim);
    S_reference = cnpy::getFromRawBuffer<double>(dim, 1, raw_S_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
    	BOOST_REQUIRE_CLOSE(S_results(i), S_reference(i), 1.0e-12);
    }
    */
    BOOST_TEST_MESSAGE("S operator diagonal by numerical quadrature");
    for (int i = 0; i < cavity.size(); ++i) {
	    BOOST_TEST_MESSAGE("S_{" << i+1 << ", " << i+1 << "} = " 
			    << std::setprecision(std::numeric_limits<long double>::digits10) << S_results(i));
    }
    
    Eigen::VectorXd D_results = Eigen::VectorXd::Zero(cavity.size());
    for (int i = 0; i < cavity.size(); ++i) {
	D_results(i) = gf.diagonalD(cavity.elements(i)); 
    }
    // In case you need to update the reference files...
    /*
    cnpy::npy_save("anisotropic_D_numerical.npy", D_results.data(), shape, 1, "w", false);
    cnpy::NpyArray raw_D_ref = cnpy::npy_load("anisotropic_D_numerical.npy");
    dim = raw_D_ref.shape[0];
    Eigen::VectorXd D_reference = Eigen::VectorXd::Zero(dim);
    D_reference = cnpy::getFromRawBuffer<double>(dim, 1, raw_D_ref.data);
    for (int i = 0; i < cavity.size(); ++i) {
    	BOOST_REQUIRE_CLOSE(D_results(i), D_reference(i), 1.0e-12);
    }
    */
    BOOST_TEST_MESSAGE("D operator diagonal by numerical quadrature");
    for (int i = 0; i < cavity.size(); ++i) {
	    BOOST_TEST_MESSAGE("D_{" << i+1 << ", " << i+1 << "} = "
			    << std::setprecision(std::numeric_limits<long double>::digits10) << D_results(i));
    }
}
