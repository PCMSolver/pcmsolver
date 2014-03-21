#ifndef TAYLORPIMPL_HPP
#define TAYLORPIMPL_HPP

// Disable obnoxious warnings from Eigen headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wextra"
#include "taylor.hpp" 
#pragma GCC diagnostic pop
#elif (__INTEL_COMPILER)
#pragma warning push
#pragma warning disable "-Wall"
#include "taylor.hpp" 
#pragma warning pop
#endif

#endif // TAYLORPIMPL_HPP
