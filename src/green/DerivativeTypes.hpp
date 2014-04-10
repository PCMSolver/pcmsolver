#ifndef DERIVATIVETYPES_HPP
#define DERIVATIVETYPES_HPP

#include "TaylorPimpl.hpp"

#include <boost/mpl/vector.hpp>

typedef taylor<double, 1, 1> AD_directional;
typedef taylor<double, 3, 1> AD_gradient;
typedef taylor<double, 3, 2> AD_hessian;

typedef boost::mpl::vector<double, AD_directional, AD_gradient, AD_hessian>
derivative_types;

#endif // DERIVATIVETYPES_HPP
