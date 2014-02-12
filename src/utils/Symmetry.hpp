#ifndef SYMMETRY_HPP
#define SYMMETRY_HPP

#include <string>
#include <vector>

#include "Config.hpp"

// Disable obnoxious warnings from Eigen headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall" 
#pragma GCC diagnostic ignored "-Weffc++" 
#pragma GCC diagnostic ignored "-Wextra"
#include <Eigen/Dense>
#pragma GCC diagnostic pop
#elif (__INTEL_COMPILER)
#pragma warning push
#pragma warning disable "-Wall"
#include <Eigen/Dense>
#pragma warning pop
#endif

/*! \file Symmetry.hpp
 *  \class Symmetry
 *  \brief Contains very basic info about symmetry (only Abelian groups)
 *  \author Roberto Di Remigio
 *  \date 2014
 *  
 *   The PCM part needs to know very little about point group symmetry:
 *   1. we need to pass the correct string with the group name to PEDRA;       
 *   2. PEDRA will construct the cavity using the symmetric generator;         
 *   3. Once we have the correct cavity we build the whole PCM matrix;         
 *   4. Last step is to block diagonalize the PCM matrix.                      
 *      The transformation is built using the concept of parity of bitstrings: 
 *       Eigen::MatrixXd U(maxrep, maxrep);                                    
 *       for ( int i = 0; i < maxrep; ++i)                                     
 *       {                                                                     
 *       		for ( int j = 0; j < maxrep; ++j)                     
 *       		{                                                     
 *       			U(i, j) = parity(i & j);                      
 *       		}                                                     
 *       }                                                                     
 *       maxrep is the number of nontrivial symmetry operations,               
 *       i.e. nr_operations - 1 (the identity - E - is the trivial             
 *       operation in every group)                                             
 *       The parity vector is defined as:                                      
 *       Eigen::VectorXd parity(8) << 1.0, -1.0, -1.0, 1.0,                    
 *                                   -1.0,  1.0, 1.0, -1.0;                    
 *                                                                             
 *   The conclusion is, instead of adapting  PSI4 and mpqc way of handling     
 *   symmetry, let us just use what we know, i.e. the DALTON (and DIRAC)       
 *   way of handling symmetry.                                                 
 *   Lots of the infrastructure will not be needed as we don't actually        
 *   need to build character tables and matrix representations.                
 */

class Symmetry
{
	private:
	public:
		static int parity(int i);
};

#endif // SYMMETRY_HPP
