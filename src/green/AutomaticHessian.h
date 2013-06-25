#ifndef AUTOMATICHESSIAN_H
#define AUTOMATICHESSIAN_H

#include "taylor.hpp"

#include "Evaluate.h"

/*! \file AutomaticHessian.h
 *  \class AutomaticHessian
 *  \brief Implements an automatic directional derivative and Hessian. 
 *  \author Roberto Di Remigio
 *  \date 2013
 *
 *  This class implements the evaluation of the function, its directional
 *  derivative and its Hessian using automatic differentiation:
 *  \f$\mathbf{\nabla}_{\mathbf{u}} f(\mathbf{x}) = \mathbf{u}^{t} \mathbf{\nabla} f(\mathbf{x})\f$, 
 *  \f$\mathbf{\nabla}_{\mathbf{u}}\mathbf{\nabla}_{\mathbf{v}}f(\mathbf{x}) = \mathbf{u}^{t} \mathbf{H}(f(\mathbf{x})) \mathbf{v}\f$
 *  The automatic differentiation engine is Ulf Ekström's libtaylor library.
 */

class AutomaticHessian : public Evaluate
{
	public:
		virtual void evaluate()
		{
		}
};

#endif // AUTOMATICHESSIAN_H
