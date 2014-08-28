#include "Vacuum.hpp"

#include <cmath>

#include "Config.hpp"

#include <Eigen/Dense>
#include "TaylorPimpl.hpp"

void Vacuum::compDiagonal(const Eigen::VectorXd & elementArea_,
                          const Eigen::VectorXd & elementRadius_, Eigen::MatrixXd & S_,
                          Eigen::MatrixXd & D_) const
{
    int size = S_.rows();

    double factor = 1.07;

    for (int i = 0; i < size; ++i) {
        S_(i, i) = factor * sqrt(4 * M_PI / elementArea_(i));
        D_(i, i) = - factor * sqrt(M_PI / elementArea_(i)) * (1 / elementRadius_(i));
    }
}

Eigen::Array4d Vacuum::numericalDirectional(Eigen::Vector3d & sourceNormal_,
        Eigen::Vector3d & source_, Eigen::Vector3d & probeNormal_,
        Eigen::Vector3d & probe_) const
{
    Eigen::Array4d result = Eigen::Array4d::Zero();

    // The finite difference step
    double delta = 1.0e-4;

    // Value of the function
    result(0) = 1 / ((source_ - probe_).norm());

    Eigen::Vector3d deltaPlus = probe_ + (probeNormal_.normalized() * delta);
    Eigen::Vector3d deltaMinus = probe_ - (probeNormal_.normalized() * delta);

    double funcPlus =  1 / ((source_ - deltaPlus).norm());
    double funcMinus = 1 / ((source_ - deltaMinus).norm());
    // Directional derivative wrt probe_
    result(1) = (funcPlus - funcMinus)/(2.0 * delta);

    deltaPlus = source_ + (sourceNormal_.normalized() * delta);
    deltaMinus = source_ - (sourceNormal_.normalized() * delta);

    funcPlus =  1 / ((deltaPlus - probe_).norm());
    funcMinus = 1 / ((deltaMinus - probe_).norm());
    // Directional derivative wrt source_
    result(2) =  (funcPlus - funcMinus)/(2.0 * delta);

    // Value of the Hessian
    result(3) = 0;

    return result;
}

Eigen::Array4d Vacuum::analyticDirectional(Eigen::Vector3d & sourceNormal_,
        Eigen::Vector3d & source_, Eigen::Vector3d & probeNormal_,
        Eigen::Vector3d & probe_) const
{
    Eigen::Array4d result = Eigen::Array4d::Zero();
    double distance = (source_ - probe_).norm();
    double distance_3 = pow(distance, 3.0);
    double distance_5 = pow(distance, 5.0);

    // Value of the function
    result(0) = 1 / distance;
    // Value of the directional derivative wrt probe_
    result(1) = (source_ - probe_).dot(probeNormal_) / distance_3 ;
    // Directional derivative wrt source_
    result(2) = - (source_ - probe_).dot(sourceNormal_) / distance_3;
    // Value of the Hessian
    result(3) = sourceNormal_.dot(probeNormal_) / distance_3 - 3 * ((
                    source_ - probe_).dot(sourceNormal_))*((source_ - probe_).dot(
                                probeNormal_)) / distance_5;

    return result;
}

Eigen::Array4d Vacuum::automaticDirectional(Eigen::Vector3d & sourceNormal_,
        Eigen::Vector3d & source_, Eigen::Vector3d & probeNormal_,
        Eigen::Vector3d & probe_) const
{
    Eigen::Array4d result = Eigen::Array4d::Zero();
    taylor<double, 3, 1> dx(0, 0), dy(0, 1), dz(0, 2);

    taylor<double, 3, 1> tmp;

    // Directional derivative wrt probe_
    tmp = 1 / sqrt((source_(0) - (probe_(0) + probeNormal_(0) * dx)) * (source_(0) -
                   (probe_(0) + probeNormal_(0) * dx))
                   + (source_(1) - (probe_(1) + probeNormal_(1) * dy)) * (source_(1) - (probe_
                           (1) + probeNormal_(1) * dy))
                   + (source_(2) - (probe_(2) + probeNormal_(2) * dz)) * (source_(2) - (probe_
                           (2) + probeNormal_(2) * dz)));
    tmp.deriv_facs();

    // Value of the function
    result(0) = tmp[0];
    // Value of the directional derivative wrt probe_
    for (int i = 1; i < 4; ++i) {
        result(1) += tmp[i];
    }

    // Directional derivative wrt source_
    tmp = 1 / sqrt(((source_(0) + sourceNormal_(0) * dx) - probe_(0)) * ((source_
                   (0) + sourceNormal_(0) * dx) - probe_(0))
                   + ((source_(1) + sourceNormal_(1) * dy) - probe_(1)) * ((source_(1) + sourceNormal_
                           (1) * dy) - probe_(1))
                   + ((source_(2) + sourceNormal_(2) * dz) - probe_(2)) * ((source_(2) + sourceNormal_
                           (2) * dz) - probe_(2)));
    tmp.deriv_facs();

    // Value of the directional derivative wrt source_
    for (int i = 1; i < 4; ++i) {
        result(2) += tmp[i];
    }
    // Value of the Hessian
    result(3) = 0;

    return result;
}


Eigen::Array4d Vacuum::automaticGradient(Eigen::Vector3d & sourceNormal_,
        Eigen::Vector3d & source_, Eigen::Vector3d & probeNormal_,
        Eigen::Vector3d & probe_) const
{
    Eigen::Array4d result = Eigen::Array4d::Zero();
    taylor<double, 3, 1> dx(0, 0), dy(0, 1), dz(0, 2);

    taylor<double, 3, 1> tmp;

    tmp = 1 / sqrt((source_(0) - (probe_(0) + dx)) * (source_(0) - (probe_(0) + dx))
                   + (source_(1) - (probe_(1) + dy)) * (source_(1) - (probe_(1) + dy))
                   + (source_(2) - (probe_(2) + dz)) * (source_(2) - (probe_(2) + dz)));
    tmp.deriv_facs();

    // Value of the function
    result(0) = tmp[0];
    // Value of the directional derivative wrt probe_
    for (int i = 1; i < 4; ++i) {
        result(1) += tmp[i] * probeNormal_(i-1);
    }
    // Value of the directional derivative wrt source_
    for (int i = 1; i < 4; ++i) {
        result(2) += -tmp[i] * sourceNormal_(i-1);
    }

    // Value of the Hessian
    result(3) = 0;

    return result;
}

Eigen::Array4d Vacuum::automaticHessian(Eigen::Vector3d & sourceNormal_,
                                        Eigen::Vector3d & source_, Eigen::Vector3d & probeNormal_,
                                        Eigen::Vector3d & probe_) const
{
    Eigen::Array4d result = Eigen::Array4d::Zero();

    taylor<double, 3, 2> dx(0, 0), dy(0, 1), dz(0, 2);

    taylor<double, 3, 2> tmp;

    //  Gradient wrt probe_
    tmp = 1 / sqrt((source_(0) - (probe_(0) + dx)) * (source_(0) - (probe_(0) + dx))
                   + (source_(1) - (probe_(1) + dy)) * (source_(1) - (probe_(1) + dy))
                   + (source_(2) - (probe_(2) + dz)) * (source_(2) - (probe_(2) + dz)));

    tmp.deriv_facs();

    // Value of the function
    result(0) = tmp[0];
    // Value of the directional derivative wrt probe_
    for (int i = 1; i < 4; ++i) {
        result(1) += tmp[i] * probeNormal_(i-1);
    }
    // Value of the directional derivative wrt source_
    for (int i = 1; i < 4; ++i) {
        result(2) += -tmp[i] * sourceNormal_(i-1);
    }

    Eigen::Matrix3d hessian = Eigen::Matrix3d::Zero();
    // Yes, this is quite clumsy...
    hessian(0, 0) = tmp[4];
    hessian(0, 1) = tmp[5];
    hessian(0, 2) = tmp[6];
    hessian(1, 1) = tmp[7];
    hessian(1, 2) = tmp[8];
    hessian(2, 2) = tmp[9];
    hessian(1, 0) = hessian(0, 1);
    hessian(2, 0) = hessian(0, 2);
    hessian(2, 1) = hessian(1, 2);

    // Value of the Hessian
    result(3) = -sourceNormal_.transpose() * hessian * probeNormal_;

    return result;
}

std::ostream & Vacuum::printGreensFunction(std::ostream & os)
{
    os << "Green's function type: vacuum";
    return os;
}
