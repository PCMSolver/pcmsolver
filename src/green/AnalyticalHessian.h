#ifndef ANALYTICALHESSIAN_H
#define ANALYTICALHESSIAN_H

#include "Evaluate.h"

/*! \file AnalyticalHessian.h
 *  \class AnalyticalHessian
 *  \brief Implements an analytical directional derivative and Hessian. 
 *  \author Roberto Di Remigio
 *  \date 2013
 *
 *  This class implements the evaluation of the function, of its directional
 *  derivative and its Hessian using analytic differentiation formulae:
 *  @f[ 
 *  \mathbf{\nabla}_{\mathbf{u}} f(\mathbf{x}) &= \mathbf{u}^{t} \mathbf{\nabla} f(\mathbf{x}) \\
 *  \mathbf{\nabla}_{\mathbf{u}}\mathbf{\nabla}_{\mathbf{v}}f(\mathbf{x}) &= \mathbf{u}^{t} \mathbf{H}(f(\mathbf{x})) \mathbf{v}
 *  @f]
 */

class AnalyticalHessian : public Evaluate
{
	public:
		virtual void evaluate()
		{
		}
};

#endif // ANALYTICALHESSIAN_H
