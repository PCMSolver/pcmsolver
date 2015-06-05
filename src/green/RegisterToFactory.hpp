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

#ifndef REGISTERTOFACTORY_HPP
#define REGISTERTOFACTORY_HPP

#include <string>

#include <Eigen/Dense>

#include "DerivativeTypes.hpp"
#include "ForIdGreen.hpp"
#include "GreenData.hpp"
#include "GreensFunctionFactory.hpp"
#include "IntegratorForward.hpp"
#include "IntegratorTypes.hpp"

/*! \file RegisterToFactory.hpp
 *  \brief Register each Green's function to the factory.
 *  \author Roberto Di Remigio
 *  \date 2015
 *
 *  This file collects all the creational functions needed for the creation
 *  of the Green's function objects by means of the factory method.
 *  Originally, each of them was in the same file as the respective class.
 *  This, however, lead to intricate inclusion dependencies.
 */

namespace
{
    struct buildVacuum
    {
        template <typename T, typename U>
        IGreensFunction * operator()(const greenData & /* data */) {
            return new Vacuum<T, U>();
        }
    };

    IGreensFunction * createVacuum(const greenData & data)
    {
        buildVacuum build;
        return for_id<derivative_types, integrator_types>(build, data, data.howDerivative, data.howIntegrator);
    }
    const std::string VACUUM("VACUUM");
    const bool registeredVacuum =
        GreensFunctionFactory::TheGreensFunctionFactory().registerGreensFunction(
            VACUUM, createVacuum);
}

namespace
{
    struct buildUniformDielectric {
        template <typename T, typename U>
        IGreensFunction * operator()(const greenData & data) {
            return new UniformDielectric<T, U>(data.epsilon);
        }
    };

    IGreensFunction * createUniformDielectric(const greenData & data)
    {
        buildUniformDielectric build;
        return for_id<derivative_types, integrator_types>(build, data, data.howDerivative, data.howIntegrator);
    }
    const std::string UNIFORMDIELECTRIC("UNIFORMDIELECTRIC");
    const bool registeredUniformDielectric =
        GreensFunctionFactory::TheGreensFunctionFactory().registerGreensFunction(
            UNIFORMDIELECTRIC, createUniformDielectric);
}

namespace
{
    struct buildIonicLiquid {
        template <typename T, typename U>
        IGreensFunction * operator()(const greenData & data) {
            return new IonicLiquid<T, U>(data.epsilon, data.kappa);
        }
    };

    IGreensFunction * createIonicLiquid(const greenData & data)
    {
        buildIonicLiquid build;
        return for_id<derivative_types, integrator_types>(build, data, data.howDerivative, data.howIntegrator);
    }
    const std::string IONICLIQUID("IONICLIQUID");
    const bool registeredIonicLiquid =
        GreensFunctionFactory::TheGreensFunctionFactory().registerGreensFunction(
            IONICLIQUID, createIonicLiquid);
}

namespace
{
    struct buildAnisotropicLiquid {
        template <typename T, typename U>
        IGreensFunction * operator()(const greenData & data) {
            return new AnisotropicLiquid<T, U>(data.epsilonTensor, data.eulerAngles);
        }
    };

    IGreensFunction * createAnisotropicLiquid(const greenData & data)
    {
        buildAnisotropicLiquid build;
        return for_id<derivative_types, integrator_types>(build, data, data.howDerivative, data.howIntegrator);
    }
    const std::string ANISOTROPICLIQUID("ANISOTROPICLIQUID");
    const bool registeredAnisotropicLiquid =
        GreensFunctionFactory::TheGreensFunctionFactory().registerGreensFunction(
            ANISOTROPICLIQUID, createAnisotropicLiquid);
}

namespace
{
    struct buildSphericalDiffuse {
        template <typename T, typename U>
        IGreensFunction * operator()(const greenData & data) {
            return new SphericalDiffuse<T, U>(data.epsilon1, data.epsilon2, data.width, data.center, data.origin);
        }
    };

    IGreensFunction * createSphericalDiffuse(const greenData & data)
    {
        buildSphericalDiffuse build;
        return for_id<integrator_types, onelayer_diffuse_profile_types>(build, data, data.howIntegrator, data.howProfile);
    }
    const std::string SPHERICALDIFFUSE("SPHERICALDIFFUSE");
    const bool registeredSphericalDiffuse =
        GreensFunctionFactory::TheGreensFunctionFactory().registerGreensFunction(
            SPHERICALDIFFUSE, createSphericalDiffuse);
}

#endif // REGISTERTOFACTORY_HPP
