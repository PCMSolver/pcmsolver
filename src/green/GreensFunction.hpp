#ifndef GREENSFUNCTION_HPP
#define GREENSFUNCTION_HPP

#include <iosfwd>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "IGreensFunction.hpp"

/*! \file GreensFunction.hpp
 *  \class GreensFunction
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
    /*! 
     *  Returns value of the Greens's function for the pair
     *  of points p1, p2: \f$ G(\mathbf{p}_1, \mathbf{p}_2)\f$
     *  
     *  \param[in] p1 first point
     *  \param[in] p2 second point
     */
    virtual double function(const Eigen::Vector3d & p1, const Eigen::Vector3d &p2) const;
    /*! 
     *  Returns value of the directional derivative of the 
     *  Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot \mathbf{n}_{\mathbf{p}_2}\f$
     *  Notice that this method returns the directional derivative with respect
     *  to the probe point, thus assuming that the direction is relative to that point.
     *  
     *  \param[in] direction the direction
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    virtual double derivative(const Eigen::Vector3d & direction,
                              const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const = 0;
    /*! 
     *  Returns value of the directional derivative of the 
     *  Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_1}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot \mathbf{n}_{\mathbf{p}_1}\f$
     *  Notice that this method returns the directional derivative with respect
     *  to the source point.
     *  
     *  \param[in] normal_p1 the normal vector to p1 
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    virtual double derivativeSource(const Eigen::Vector3d & normal_p1,
                                    const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const;
    /*! 
     *  Returns value of the directional derivative of the 
     *  Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\cdot \mathbf{n}_{\mathbf{p}_2}\f$
     *  Notice that this method returns the directional derivative with respect
     *  to the probe point.
     *  
     *  \param[in] normal_p2 the normal vector to p1 
     *  \param[in]        p1 first point
     *  \param[in]        p2 second point
     */
    virtual double derivativeProbe(const Eigen::Vector3d & normal_p2,
                                   const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const;
    /*! 
     *  Returns full gradient of Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_1}}G(\mathbf{p}_1, \mathbf{p}_2)\f$
     *  Notice that this method returns the gradient with respect to the source point.
     *  
     *  \param[in] p1 first point
     *  \param[in] p2 second point
     */
    virtual Eigen::Vector3d gradientSource(const Eigen::Vector3d & p1,
                                           const Eigen::Vector3d & p2) const;
    /*! 
     *  Returns full gradient of Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\f$
     *  Notice that this method returns the gradient with respect to the probe point.
     *  
     *  \param[in] p1 first point
     *  \param[in] p2 second point
     */
    virtual Eigen::Vector3d gradientProbe(const Eigen::Vector3d & p1,
                                          const Eigen::Vector3d & p2) const;
    /*! 
     *  Returns full gradient of Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_1}}G(\mathbf{p}_1, \mathbf{p}_2)\f$
     *  Notice that this method returns the gradient with respect to the source point.
     * 
     *  \param[in] gradient the gradient
     *  \param[in]       p1 first point
     *  \param[in]       p2 second point
     */
    virtual void gradientSource(Eigen::Vector3d & gradient, const Eigen::Vector3d & p1,
                                const Eigen::Vector3d & p2) const;
    /*! 
     *  Returns full gradient of Greens's function for the pair of points p1, p2:
     *  \f$ \nabla_{\mathbf{p_2}}G(\mathbf{p}_1, \mathbf{p}_2)\f$
     *  Notice that this method returns the gradient with respect to the probe point.
     *
     *  \param[in] gradient the gradient
     *  \param[in]       p1 first point
     *  \param[in]       p2 second point
     */
    virtual void gradientProbe(Eigen::Vector3d & gradient, const Eigen::Vector3d & p1,
                               const Eigen::Vector3d & p2) const;
    
    virtual void delta(double value);
    virtual double delta() { return delta_; }
    virtual double epsilon() const = 0;
    
    friend std::ostream & operator<<(std::ostream & os, GreensFunction & gf) {
        return gf.printObject(os);
    }
protected:
    /*!
     *  Evaluates the Green's function given a pair of points
     *
     *  \param[in] source the source point
     *  \param[in]  probe the probe point
     */
    virtual T operator()(T * source, T * probe) const = 0;
    virtual std::ostream & printObject(std::ostream & os);
    double delta_;
};

#endif // GREENSFUNCTION_HPP
