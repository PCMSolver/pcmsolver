/*! \file PCMSolver.cpp 
\brief PCM solver
*/

#include <string>
#include <vector>
#include <iostream>
#include "PCMSolver.h"
#include "GreensFunction.h"
#include "UniformDielectric.h"
#include "GreensFunctionSum.h"
#include "MetalSphere.h"

template <class Green>
Green& PCMSolver<Green>::getGreensFunction(){
	return *greensFunction;
}

template class PCMSolver<UniformDielectric>;
template class PCMSolver<GreensFunctionSum>;
template class PCMSolver<MetalSphere>;


