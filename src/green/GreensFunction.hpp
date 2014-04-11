#ifndef GREENSFUNCTION_HPP
#define GREENSFUNCTION_HPP

#include <iosfwd>

#include "EigenPimpl.hpp"

#include "IGreensFunction.hpp"

/*! \file GreensFunction.hpp
 *  \class GreensFunction<T>
 *  \brief Templated interface for Green's functions
 *  \author Luca Frediani and Roberto Di Remigio
 *  \date 2012-2014
 *  \tparam T evaluation strategy for the function and its derivatives
 */

template <typename T>
class GreensFunction: public IGreensFunction
{
public:
    GreensFunction() : IGreensFunction(true), delta_(1.0e-4) {}
    GreensFunction(bool uniform) : IGreensFunction(uniform), delta_(1.0e-4) {}
    virtual ~GreensFunction() {}
    virtual double function(const Eigen::Vector3d & p1, const Eigen::Vector3d &p2) const;
    virtual double derivative(const Eigen::Vector3d & direction,
                              const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const = 0;
    virtual double derivativeSource(const Eigen::Vector3d & direction,
                                    const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const;
    virtual double derivativeProbe(const Eigen::Vector3d & direction,
                                   const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const;
    virtual Eigen::Vector3d gradientSource(const Eigen::Vector3d & p1,
                                           const Eigen::Vector3d & p2) const;
    virtual Eigen::Vector3d gradientProbe(const Eigen::Vector3d & p1,
                                          const Eigen::Vector3d & p2) const;
    virtual void gradientSource(Eigen::Vector3d & gradient, const Eigen::Vector3d & p1,
                                const Eigen::Vector3d & p2) const;
    virtual void gradientProbe(Eigen::Vector3d & gradient, const Eigen::Vector3d & p1,
                               const Eigen::Vector3d & p2) const;
    virtual void delta(double value);
    virtual double delta() { return delta_; }
    virtual double epsilon() const = 0;
    virtual double compDiagonalElementS(double area) const = 0;
    virtual double compDiagonalElementD(double area, double radius) const = 0;
    /*!
     *  Get the S and D matrices
     */
    virtual void operator()(Eigen::MatrixXd & S, Eigen::MatrixXd & D,
                            const Eigen::MatrixXd & centers, const Eigen::MatrixXd & normals,
                            const Eigen::VectorXd & areas, const Eigen::VectorXd & radii) const = 0;
    /*!
     *  Get the S matrix
     */
    virtual void operator()(Eigen::MatrixXd & S,
                            const Eigen::MatrixXd & centers, const Eigen::MatrixXd & normals,
                            const Eigen::VectorXd & areas) const = 0;
    /*!
     *  Get the D matrix
     */
    virtual void operator()(Eigen::MatrixXd & D,
                            const Eigen::MatrixXd & centers, const Eigen::MatrixXd & normals,
                            const Eigen::VectorXd & areas, const Eigen::VectorXd & radii) const = 0;

    friend std::ostream & operator<<(std::ostream & os, GreensFunction & gf) {
        return gf.printObject(os);
    }
protected:
    virtual T evaluate(T * source, T * probe) const = 0;
    virtual std::ostream & printObject(std::ostream & os);
    double delta_;
};

#endif // GREENSFUNCTION_HPP
