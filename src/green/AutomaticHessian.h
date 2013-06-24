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
 *  derivative and its Hessian using automatic differentiation.
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
