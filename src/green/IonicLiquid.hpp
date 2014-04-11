#ifndef IONICLIQUID_HPP
#define IONICLIQUID_HPP

#include <iosfwd>
#include <string>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "DerivativeTypes.hpp"
#include "ForIdGreen.hpp"
#include "GreenData.hpp"
#include "GreensFunction.hpp"
#include "GreensFunctionFactory.hpp"
#include "IGreensFunction.hpp"

/*! \file IonicLiquid.hpp
 *  \class IonicLiquid<T>
 *  \brief Green's functions for ionic liquid, described by the linearized Poisson-Boltzmann equation.
 *  \author Luca Frediani, Roberto Di Remigio
 *  \date 2013-2014
 *  \tparam T evaluation strategy for the function and its derivatives
 */

template <typename T>
class IonicLiquid : public GreensFunction<T>
{
public:
    IonicLiquid(double epsilon, double kappa) : GreensFunction<T>(false),
        epsilon_(epsilon), kappa_(kappa) {}
    virtual ~IonicLiquid() {}
    virtual double derivative(const Eigen::Vector3d & direction,
                              const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const;
    virtual double compDiagonalElementS(double area) const ;
    virtual double compDiagonalElementD(double area, double radius) const;
    virtual double epsilon() const { return epsilon_; }
    /*!
     *  Get the S and D matrices
     */
    virtual void operator()(Eigen::MatrixXd & S, Eigen::MatrixXd & D,
                            const Eigen::MatrixXd & centers, const Eigen::MatrixXd & normals,
                            const Eigen::VectorXd & areas, const Eigen::VectorXd & radii) const;
    /*!
     *  Get the S matrix
     */
    virtual void operator()(Eigen::MatrixXd & S,
                            const Eigen::MatrixXd & centers, const Eigen::MatrixXd & normals,
                            const Eigen::VectorXd & areas) const;
    /*!
     *  Get the D matrix
     */
    virtual void operator()(Eigen::MatrixXd & D,
                            const Eigen::MatrixXd & centers, const Eigen::MatrixXd & normals,
                            const Eigen::VectorXd & areas, const Eigen::VectorXd & radii) const;

    friend std::ostream & operator<<(std::ostream & os, IonicLiquid & gf) {
        return gf.printObject(os);
    }
private:
    virtual T evaluate(T * source, T * probe) const;
    double epsilon_;
    double kappa_;
    virtual std::ostream & printObject(std::ostream & os);
};

namespace
{
    struct buildIonicLiquid {
        template <typename DerivativeType>
        IGreensFunction * operator()(const greenData & _data) {
            return new IonicLiquid<DerivativeType>(_data.epsilon, _data.kappa);
        }
    };

    IGreensFunction * createIonicLiquid(const greenData & _data)
    {
        buildIonicLiquid build;
        return for_id<derivative_types>(build, _data, _data.how);
    }
    const std::string IONICLIQUID("IonicLiquid");
    const bool registeredIonicLiquid =
        GreensFunctionFactory::TheGreensFunctionFactory().registerGreensFunction(
            IONICLIQUID, createIonicLiquid);
}

#endif // IONICLIQUID_HPP
