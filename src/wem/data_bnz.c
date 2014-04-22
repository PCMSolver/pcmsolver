/* warning-disabler-start */

#if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wextra"
#elif defined(__ICC) || defined(__INTEL_COMPILER)
#pragma warning push
#pragma warning disable "-Wall"
#elif defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wall"
#pragma clang diagnostic ignored "-Weffc++"
#pragma clang diagnostic ignored "-Wextra"
#pragma clang diagnostic ignored "-Wdocumentation"
#pragma clang diagnostic ignored "-Wdeprecated-register"
#pragma clang diagnostic ignored "-Wincompatible-pointer-types"
#pragma clang diagnostic ignored "-Wempty-body"
#endif

/* warning-disabler-end */

const unsigned int k = 12;
const vector3 x[12] = {
    {5.274, 1.999, -8.568},
    {6.627, 2.018, -8.209},
    {7.366, 0.829, -8.202},
    {6.752, -0.379, -8.554},
    {5.399, -0.398, -8.912},
    {4.660, 0.791, -8.919},
    {4.704, 2.916, -8.573},
    {7.101, 2.950, -7.938},
    {8.410, 0.844, -7.926},
    {7.322, -1.296, -8.548},
    {4.925, -1.330, -9.183},
    {3.616, 0.776, -9.196}
};
const double alpha[12] = { 6, 6, 6, -6, -6, -6, 1, 1, 1, -1, -1, 1 };
/* warning-disabler-start */

#if (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
#pragma GCC diagnostic pop
#elif defined(__ICC) || defined(__INTEL_COMPILER)
#pragma warning pop
#elif defined(__clang__)
#pragma clang diagnostic pop
#endif

/* warning-disabler-end */

