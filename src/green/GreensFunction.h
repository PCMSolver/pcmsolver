#ifndef GREENSFUNCTION_H
#define GREENSFUNCTION_H

/*! \file GreensFunction.h
 *  \class GreensFunction
 *  \brief Abstract base class for the Green's function generator. 
 *  \author Luca Frediani
 *  \date 2011
 *  
 *  A generic Green´s function to represent the electrostatic potential for a given environment
 */

#include "GreensFunctionInterface.h"

class Section;

template<typename T>
class GreensFunction: public GreensFunctionInterface
{
 public:
    GreensFunction(){delta = 1.0e-4;}
    virtual ~GreensFunction(){};

    // From GreensFunctionInterface
    /*! 
     * \brief Value of the Green's function for a pair of source and probe points.
     * \param[in] p1 the source point.
     * \param[in] p2 the probe point.
     *
     * This function is used in computing the kernel of the \f$\mathcal{S}\f$ integral operator.
     */
    virtual double evalf(Eigen::Vector3d &p1, Eigen::Vector3d &p2);
    /*!
     * \brief Wrapper function for the evaluation of the derivative of the Green's function. 
     * \param[in] direction the direction used to calculate the directional derivative.
     * \param[in] p1 the source point.
     * \param[in] p2 the probe point.
     *
     * This function is used in computing the kernel of both the \f$\mathcal{D}\f$ and 
     * \f$\mathcal{D}^\dagger\f$ integral operators. It is a wrapper for the other derivative calculation
     * functions provided in this Abstract Base Class.
     */
    virtual double evald(Eigen::Vector3d &direction, Eigen::Vector3d &p1, Eigen::Vector3d &p2) = 0;
    /*!
     * \brief Directional derivative of the Green's function in a direction relative to the source point.
     * \param[in] direction the direction used to calculate the directional derivative.
     * \param[in] p1 the source point.
     * \param[in] p2 the probe point. 
     *
     * This function is used in computing the kernel of the \f$\mathcal{D}^\dagger\f$ integral operator.
     */
    virtual double derivativeSource(Eigen::Vector3d &direction, Eigen::Vector3d &p1, Eigen::Vector3d &p2);
    /*!
     * \brief Directional derivative of the Green's function in a direction relative to the probe point.
     * \param[in] direction the direction used to calculate the directional derivative.
     * \param[in] p1 the source point.
     * \param[in] p2 the probe point.
     *
     * This function is used in computing the kernel of the \f$\mathcal{D}\f$ integral operator.
     */
    virtual double derivativeProbe(Eigen::Vector3d &direction, Eigen::Vector3d &p1, Eigen::Vector3d &p2);

    /*!
     * \brief Gradient of the Green's function with respect to the source point.
     * \param[in] p1 the source point.
     * \param[in] p2 the probe point. 
     *
     * This function is used in computing the kernel of the \f$\mathcal{D}^\dagger\f$ integral operator.
     */
    virtual Eigen::Vector3d gradientSource(Eigen::Vector3d &p1, Eigen::Vector3d &p2);
    /*!
     * \brief Gradient of the Green's function with respect to the probe point.
     * \param[in] p1 the source point.
     * \param[in] p2 the probe point.
     *
     * This function is used in computing the kernel of the \f$\mathcal{D}\f$ integral operator.
     */
    virtual Eigen::Vector3d gradientProbe(Eigen::Vector3d &p1, Eigen::Vector3d &p2);
    virtual double getDielectricConstant();
    /*!
     * \brief Gradient of the Green's function with respect to the source point.
     * \param[out] gradient a vector containing the gradient.
     * \param[in] p1 the source point.
     * \param[in] p2 the probe point.
     *
     * This function is used in computing the kernel of the \f$\mathcal{D}^\dagger\f$ integral operator.
     * It is an overloaded version using the pass-by-reference semantics.
     */
    virtual void gradientSource(Eigen::Vector3d &gradient, Eigen::Vector3d &p1, Eigen::Vector3d &p2);
    /*!
     * \brief Gradient of the Green's function with respect to the probe point.
     * \param[out] gradient a vector containing the gradient.
     * \param[in] p1 the source point.
     * \param[in] p2 the probe point.
     *
     * This function is used in computing the kernel of the \f$\mathcal{D}\f$ integral operator.
     * It is an overloaded version using the pass-by-reference semantics.
     */
    virtual void gradientProbe(Eigen::Vector3d &gradient, Eigen::Vector3d &p1, Eigen::Vector3d &p2);
    void setDelta(double value);
    double getDelta(){return delta;}
    GreensFunction<T> * allocateGreensFunction(const Section &green);
    GreensFunction<T> * allocateGreensFunction(double dielConst);
    GreensFunction<T> * allocateGreensFunction();
    virtual double compDiagonalElementS(double area) = 0;
    virtual double compDiagonalElementD(double area, double radius) = 0;
 protected:
    virtual T evalGreensFunction(T * source, T * probe) = 0;
    std::ostream & printObject(std::ostream & os);
    double delta;
    bool uniformFlag;
};

#endif // GREENSFUNCTION_H
