#ifndef GETKWPIMPL_HPP
#define GETKWPIMPL_HPP

// Disable obnoxious warnings from Getkw headers
#if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wextra"
#include "Getkw.h"
#pragma GCC diagnostic pop
#elif defined(__ICC) || defined(__INTEL_COMPILER)
#pragma warning push
#pragma warning disable "-Wall"
#include "Getkw.h"
#pragma warning pop
#elif defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wall"
#pragma clang diagnostic ignored "-Weffc++"
#pragma clang diagnostic ignored "-Wextra"
#pragma clang diagnostic ignored "-Wdocumentation"
#pragma clang diagnostic ignored "-Wdeprecated-register"
#include "Getkw.h"
#pragma clang diagnostic pop
#endif

#endif // GETKWPIMPL_HPP
