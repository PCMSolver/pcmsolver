#ifndef AUTOMATICDIRECTIONAL_H
#define AUTOMATICDIRECTIONAL_H

#include "taylor.hpp"

#include "Evaluate.h"

/*! \file AutomaticDirectional.h
 *  \class AutomaticDirectional
 *  \brief Implements an automatic directional derivative. 
 *  \author Roberto Di Remigio
 *  \date 2013
 *
 *  This class implements the evaluation of the function and of its directional
 *  derivative using automatic differentiation:
 *  \f$ \mathbf{\nabla}_{\mathbf{u}} f(\mathbf{x}) = \mathbf{u}^{t}  \mathbf{\nabla} f(\mathbf{x}) \f$
 *  The automatic differentiation engine is Ulf Ekström's libtaylor library.
 */

class AutomaticDirectional : public Evaluate
{
	public:
		virtual void evaluate()
		{
		}
};

#endif // AUTOMATICDIRECTIONAL_H
