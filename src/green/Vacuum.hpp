#ifndef VACUUM_HPP
#define VACUUM_HPP

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

/*! \file Vacuum.hpp
 *  \class Vacuum<T>
 *  \brief Green's function for vacuum.
 *  \author Luca Frediani and Roberto Di Remigio
 *  \date 2012-2014
 *  \tparam T evaluation strategy for the function and its derivatives
 */

template <typename T>
class Vacuum : public GreensFunction<T>
{
public:
    Vacuum() : GreensFunction<T>(), epsilon_(1.0) {}
    virtual ~Vacuum() {}
    virtual double derivative(const Eigen::Vector3d & direction,
                              const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const;
    virtual double compDiagonalElementS(double area) const ;
    virtual double compDiagonalElementD(double area, double radius) const;
    virtual double epsilon() const { return 1.0; }
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

    friend std::ostream & operator<<(std::ostream & os, Vacuum & gf) {
        return gf.printObject(os);
    }
private:
    virtual T evaluate(T * source, T * probe) const;
    virtual std::ostream & printObject(std::ostream & os);
    double epsilon_;
};

namespace
{
    struct buildVacuum {
        template <typename DerivativeType>
        IGreensFunction * operator()(const greenData & _data) {
            return new Vacuum<DerivativeType>();
        }
    };

    IGreensFunction * createVacuum(const greenData & _data)
    {
        buildVacuum build;
        return for_id<derivative_types>(build, _data, _data.how);
    }
    const std::string VACUUM("Vacuum");
    const bool registeredVacuum =
        GreensFunctionFactory::TheGreensFunctionFactory().registerGreensFunction(
            VACUUM, createVacuum);
}

#endif // VACUUM_HPP
