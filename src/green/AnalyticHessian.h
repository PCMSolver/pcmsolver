#ifndef ANALYTICHESSIAN_H
#define ANALYTICHESSIAN_H

#include "Evaluate.h"

/*! \file AnalyticHessian.h
 *  \class AnalyticHessian
 *  \brief Implements an analytic directional derivative and Hessian. 
 *  \author Roberto Di Remigio
 *  \date 2013
 *
 *  This class implements the evaluation of the function, of its directional
 *  derivative and its Hessian using analytic differentiation formulae:
 *  \f$\mathbf{\nabla}_{\mathbf{u}} f(\mathbf{x}) = \mathbf{u}^{t} \mathbf{\nabla} f(\mathbf{x})\f$, 
 *  \f$\mathbf{\nabla}_{\mathbf{u}}\mathbf{\nabla}_{\mathbf{v}}f(\mathbf{x}) = \mathbf{u}^{t} \mathbf{H}(f(\mathbf{x})) \mathbf{v}\f$
 */

class AnalyticHessian : public Evaluate
{
	public:
		virtual void evaluate()
		{
		}
};

#endif // ANALYTICHESSIAN_H
