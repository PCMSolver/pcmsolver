#ifndef NUMERICALINTEGRATOR_HPP
#define NUMERICALINTEGRATOR_HPP

#include <algorithm>
#include <cmath>
#include <functional>
#include <iosfwd>
#include <stdexcept>

#include "Config.hpp"

#include <Eigen/Dense>

#include "Element.hpp"
#include "MathUtils.hpp"
#include "QuadratureRules.hpp"
#include "Sphere.hpp"
#include "DerivativeTypes.hpp"
#include "AnisotropicLiquid.hpp"
#include "IonicLiquid.hpp"
#include "UniformDielectric.hpp"
#include "Vacuum.hpp"

typedef std::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)>
singleLayerIntegrand;

typedef std::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &,
		        const Eigen::Vector3d &)> doubleLayerIntegrand;

/*! \file NumericalIntegrator.hpp
 *  \class NumericalIntegrator
 *  \brief Implementation of diagonal elements of S and D using numerical integration
 *  \author Roberto Di Remigio
 *  \date 2014
 *
 *  Calculates the diagonal elements of S and D by collocation, using numerical
 *  integration.
 */

template <typename DerivativeTraits,
          typename ProfilePolicy>
struct NumericalIntegrator
{
    double computeS(const Vacuum<DerivativeTraits, NumericalIntegrator> & gf, const Element & e) const {
        using namespace std::placeholders;
	    singleLayerIntegrand F = std::bind(&Vacuum<DerivativeTraits, NumericalIntegrator>::kernelS, gf, _1, _2);
        return integrator<32, 16>(F, e);
    }
    double computeD(const Vacuum<DerivativeTraits, NumericalIntegrator> & gf, const Element & e) const {
        using namespace std::placeholders;
	    doubleLayerIntegrand F = std::bind(&Vacuum<DerivativeTraits, NumericalIntegrator>::kernelD, gf, _1, _2, _3);
        return integrator<32, 16>(F, e);
    }

    double computeS(const UniformDielectric<DerivativeTraits, NumericalIntegrator> & gf, const Element & e) const {
        using namespace std::placeholders;
	    singleLayerIntegrand F = std::bind(&UniformDielectric<DerivativeTraits, NumericalIntegrator>::kernelS, gf, _1, _2);
        return integrator<32, 16>(F, e);
    }
    double computeD(const UniformDielectric<DerivativeTraits, NumericalIntegrator> & gf, const Element & e) const {
        using namespace std::placeholders;
	    doubleLayerIntegrand F = std::bind(&UniformDielectric<DerivativeTraits, NumericalIntegrator>::kernelD, gf, _1, _2, _3);
        return integrator<32, 16>(F, e);
    }

    double computeS(const IonicLiquid<DerivativeTraits, NumericalIntegrator> & gf, const Element & e) const {
        using namespace std::placeholders;
	    singleLayerIntegrand F = std::bind(&IonicLiquid<DerivativeTraits, NumericalIntegrator>::kernelS, gf, _1, _2);
        return integrator<32, 16>(F, e);
    }
    double computeD(const IonicLiquid<DerivativeTraits, NumericalIntegrator> & gf, const Element & e) const {
        using namespace std::placeholders;
	    doubleLayerIntegrand F = std::bind(&IonicLiquid<DerivativeTraits, NumericalIntegrator>::kernelD, gf, _1, _2, _3);
        return integrator<32, 16>(F, e);
    }

    double computeS(const AnisotropicLiquid<DerivativeTraits, NumericalIntegrator> & gf, const Element & e) const {
        using namespace std::placeholders;
	    singleLayerIntegrand F = std::bind(&AnisotropicLiquid<DerivativeTraits, NumericalIntegrator>::kernelS, gf, _1, _2);
        return integrator<32, 16>(F, e);
    }
    double computeD(const AnisotropicLiquid<DerivativeTraits, NumericalIntegrator> & gf, const Element & e) const {
        using namespace std::placeholders;
	    doubleLayerIntegrand F = std::bind(&AnisotropicLiquid<DerivativeTraits, NumericalIntegrator>::kernelD, gf, _1, _2, _3);
        return integrator<32, 16>(F, e);
    }
};

#endif // NUMERICALINTEGRATOR_HPP
