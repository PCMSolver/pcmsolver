#ifndef GETKWPIMPL_HPP
#define GETKWPIMPL_HPP

// Disable obnoxious warnings from Getkw headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wextra"
#include "Getkw.h"
#pragma GCC diagnostic pop
#elif (__INTEL_COMPILER)
#pragma warning push
#pragma warning disable "-Wall"
#include "Getkw.h"
#pragma warning pop
#endif

#endif // GETKWPIMPL_HPP
