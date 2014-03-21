#ifndef CNPYPIMPL_HPP
#define CNPYPIMPL_HPP

// Disable obnoxious warnings from cnpy header
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wextra"
#include "cnpy.hpp" 
#pragma GCC diagnostic pop
#elif (__INTEL_COMPILER)
#pragma warning push
#pragma warning disable "-Wall"
#include "cnpy.hpp" 
#pragma warning pop
#endif

#endif // CNPYPIMPL_HPP
