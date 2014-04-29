#include <cmath>
#include <iostream>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include <boost/numeric/quadrature/kronrodgauss.hpp>
#include <boost/numeric/quadrature/error_estimator.hpp>
#include <boost/lambda/lambda.hpp>

#include "gtestPimpl.hpp"

namespace quadrature = boost::numeric::quadrature;

/*! \test \b NumericalQuadrature_x2sinx tests adaptive integration of x^2*sin(x) in [0, 5] 
 */
TEST(NumericalQuadrature, x2sinx)
{
  struct f
  {
    double operator()(double x) const { return (std::pow(x, 2) * std::sin(x)); }
  };
  double upper = 5.0;
  double lower = 0.0;
  double expected_result = -2.0 + 10.0 * std::sin(upper) - 23.0 * std::cos(upper);
  double actual_result;

  // declare the kernel
  quadrature::kronrod_gauss<15> integrator;

  // integrate x^2*sin(x) on [0,5]
  integrator(f(), lower, upper, actual_result);

  // display the result
  std::cout << "integral(x^2*sin(x)) on [0, 5] is " << actual_result << std::endl;
  EXPECT_DOUBLE_EQ(expected_result, actual_result);
}
