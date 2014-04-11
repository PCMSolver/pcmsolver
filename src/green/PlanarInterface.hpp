#ifndef PLANARINTERFACE_HPP
#define PLANARINTERFACE_HPP

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

/*! \file PlanarInterface.hpp
 *  \class PlanarInterface<T>
 *  \brief Green's functions for a planar interface 
 *  \author Luca Frediani, Roberto Di Remigio
 *  \date 2013-2014
 *  \tparam T evaluation strategy for the function and its derivatives
 *  
 *  Reference:
 *  http://dx.doi.org/10.1063/1.1643727 
 */

template <typename T>
class PlanarInterface : public GreensFunction<T>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // See http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    PlanarInterface(double eps1, double eps2, const Eigen::Vector3d & pos, double width)
        : GreensFunction<T>(false), eps1_(eps1), eps2_(eps2), pos_(pos), width_(width), computed_(false) {}
    virtual ~PlanarInterface() {}
    virtual double derivative(const Eigen::Vector3d & direction,
                              const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const;
    virtual double compDiagonalElementS(double area) const ;
    virtual double compDiagonalElementD(double area, double radius) const;
    virtual double epsilon() const { return eps1_; } // This is just to get it to compile...
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

    friend std::ostream & operator<<(std::ostream & os, PlanarInterface & gf) {
        return gf.printObject(os);
    }
private:
    virtual T evaluate(T * source, T * probe) const;
    double eps1_;
    double eps2_;
    Eigen::Vector3d pos_;
    double width_;
    bool computed_;
    virtual std::ostream & printObject(std::ostream & os);
};

namespace
{
    struct buildPlanarInterface {
        template <typename DerivativeType>
        IGreensFunction * operator()(const greenData & _data) {
	    // We pass some bogus arguments...
	    Eigen::Vector3d orig;
	    orig << 0.0, 0.0, 0.0;
            return new PlanarInterface<DerivativeType>(_data.epsilon, 0.0, orig, 0.0);
        }
    };

    IGreensFunction * createPlanarInterface(const greenData & _data)
    {
        buildPlanarInterface build;
        return for_id<derivative_types>(build, _data, _data.how);
    }
    const std::string PLANARINTERFACE("PlanarInterface");
    const bool registeredPlanarInterface =
        GreensFunctionFactory::TheGreensFunctionFactory().registerGreensFunction(
            PLANARINTERFACE, createPlanarInterface);
}

#endif // PLANARINTERFACE_HPP
