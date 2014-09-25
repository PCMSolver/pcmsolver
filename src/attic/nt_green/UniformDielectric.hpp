#ifndef UNIFORMDIELECTRIC_HPP
#define UNIFORMDIELECTRIC_HPP

#include <iosfwd>
#include <string>

#include "Config.hpp"

#include <Eigen/Dense>

#include "GreenData.hpp"
#include "GreensFunction.hpp"
#include "GreensFunctionFactory.hpp"

/*! \file UniformDielectric.hpp
 *  \class UniformDielectric
 *  \brief Green's functions for uniform dielectric.
 *  \author Luca Frediani, Roberto Di Remigio
 *  \date 2013
 */

class UniformDielectric : public GreensFunction
{
public:
    UniformDielectric(int how_, double epsilon_) : GreensFunction(how_, true),
        epsilon(epsilon_) {}
    virtual ~UniformDielectric() {}
    virtual void compDiagonal(const Eigen::VectorXd & elementArea_,
                              const Eigen::VectorXd & elementRadius_,
                              Eigen::MatrixXd & S_, Eigen::MatrixXd & D_) const;
    virtual double dielectricConstant() const {
        return epsilon;
    }
    friend std::ostream & operator<<(std::ostream & os, UniformDielectric & green) {
        return green.printGreensFunction(os);
    }
private:
    double epsilon;
    virtual Eigen::Array4d numericalDirectional(Eigen::Vector3d & sourceNormal_,
            Eigen::Vector3d & source_,
            Eigen::Vector3d & probeNormal_, Eigen::Vector3d & probe_) const;
    virtual Eigen::Array4d analyticDirectional(Eigen::Vector3d & sourceNormal_,
            Eigen::Vector3d & source_,
            Eigen::Vector3d & probeNormal_, Eigen::Vector3d & probe_) const;
    virtual Eigen::Array4d automaticDirectional(Eigen::Vector3d & sourceNormal_,
            Eigen::Vector3d & source_,
            Eigen::Vector3d & probeNormal_, Eigen::Vector3d & probe_) const;
    virtual Eigen::Array4d automaticGradient(Eigen::Vector3d & sourceNormal_,
            Eigen::Vector3d & source_,
            Eigen::Vector3d & probeNormal_, Eigen::Vector3d & probe_) const;
    virtual Eigen::Array4d automaticHessian(Eigen::Vector3d & sourceNormal_,
                                            Eigen::Vector3d & source_,
                                            Eigen::Vector3d & probeNormal_, Eigen::Vector3d & probe_) const;
    virtual std::ostream & printGreensFunction(std::ostream & os);
};

namespace
{
    GreensFunction * createUniformDielectric(const greenData & _data)
    {
        return new UniformDielectric(_data.how, _data.epsilon);
    }
    const std::string UNIFORMDIELECTRIC("UniformDielectric");
    const bool registeredUniformDielectric =
        GreensFunctionFactory::TheGreensFunctionFactory().registerGreensFunction(
            UNIFORMDIELECTRIC, createUniformDielectric);
}

#endif // UNIFORMDIELECTRIC_HPP
