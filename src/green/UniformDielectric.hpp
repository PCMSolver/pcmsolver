#ifndef UNIFORMDIELECTRIC_HPP
#define UNIFORMDIELECTRIC_HPP

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

/*! \file UniformDielectric.hpp
 *  \class UniformDielectric<T>
 *  \brief Green's function for uniform dielectric.
 *  \author Luca Frediani and Roberto Di Remigio
 *  \date 2012-2014
 *  \tparam T evaluation strategy for the function and its derivatives
 */

template <typename T>
class UniformDielectric : public GreensFunction<T>
{
public:
    explicit UniformDielectric(double eps) : GreensFunction<T>(), epsilon_(eps) {}
    virtual ~UniformDielectric() {}
    virtual double derivative(const Eigen::Vector3d & direction,
                              const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const;
    virtual void epsilon(double eps) { epsilon_ = eps; }
    virtual double epsilon() const { return epsilon_; }
    virtual double compDiagonalElementS(double area) const ;
    virtual double compDiagonalElementD(double area, double radius) const;
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

    friend std::ostream & operator<<(std::ostream & os, UniformDielectric & gf) {
        return gf.printObject(os);
    }
private:
    virtual T evaluate(T * source, T * probe) const;
    virtual std::ostream & printObject(std::ostream & os);
    double epsilon_;
};

namespace
{
    struct buildUniformDielectric {
        template <typename DerivativeType>
        IGreensFunction * operator()(const greenData & _data) {
            return new UniformDielectric<DerivativeType>(_data.epsilon);
        }
    };

    IGreensFunction * createUniformDielectric(const greenData & _data)
    {
        buildUniformDielectric build;
        return for_id<derivative_types>(build, _data, _data.how);
    }
    const std::string UNIFORMDIELECTRIC("UniformDielectric");
    const bool registeredUniformDielectric =
        GreensFunctionFactory::TheGreensFunctionFactory().registerGreensFunction(
            UNIFORMDIELECTRIC, createUniformDielectric);
}

#endif // UNIFORMDIELECTRIC_HPP
