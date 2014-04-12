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
                              const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const;

    virtual double epsilon() const { return 1.0; }

    friend std::ostream & operator<<(std::ostream & os, Vacuum & gf) {
        return gf.printObject(os);
    }
private:
    /*!
     *  Evaluates the Green's function given a pair of points
     *
     *  \param[in] source the source point
     *  \param[in]  probe the probe point
     */
    virtual T operator()(T * source, T * probe) const;
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
