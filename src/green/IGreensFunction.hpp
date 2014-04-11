#ifndef IGREENSFUNCTION_HPP
#define IGREENSFUNCTION_HPP

#include <iosfwd>

#include "EigenPimpl.hpp"

/*! \file IGreensFunction.hpp
 *  \class IGreensFunction
 *  \brief Interface for Green's function classes
 *  \author Luca Frediani and Roberto Di Remigio
 *  \date 2012-2014
 */

class IGreensFunction
{
public:
    IGreensFunction() : uniform_(false) {}
    IGreensFunction(bool uniform) : uniform_(uniform) {}
    virtual ~IGreensFunction() {}
    virtual double function(const Eigen::Vector3d & p1,
                            const Eigen::Vector3d &p2) const = 0;
    virtual double derivative(const Eigen::Vector3d & direction,
                              const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const = 0;
    virtual double derivativeSource(const Eigen::Vector3d & direction,
                                    const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const = 0;
    virtual double derivativeProbe(const Eigen::Vector3d & direction,
                                   const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const = 0;
    virtual Eigen::Vector3d gradientSource(const Eigen::Vector3d & p1,
                                           const Eigen::Vector3d & p2) const = 0;
    virtual Eigen::Vector3d gradientProbe(const Eigen::Vector3d & p1,
                                          const Eigen::Vector3d & p2) const = 0;
    virtual void gradientSource(Eigen::Vector3d & gradient, const Eigen::Vector3d & p1,
                                const Eigen::Vector3d & p2) const = 0;
    virtual void gradientProbe(Eigen::Vector3d & gradient, const Eigen::Vector3d & p1,
                               const Eigen::Vector3d & p2) const = 0;
    virtual double epsilon() const = 0;
    virtual double compDiagonalElementS(double area) const = 0;
    virtual double compDiagonalElementD(double area, double radius) const = 0;
    bool isUniform() const { return uniform_; }
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

    friend std::ostream & operator<<(std::ostream & os, IGreensFunction & gf) {
        return gf.printObject(os);
    }
protected:
    virtual std::ostream & printObject(std::ostream & os) = 0;
    bool uniform_;
};

#endif // IGREENSFUNCTION_HPP
