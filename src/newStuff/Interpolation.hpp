#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP
/**
 * @file Interpolation.hpp
 *
 * @brief class for calculating the surface interpolation
 */
#include "Vector3.hpp"
#include "Vector2.hpp"
#include <stdlib.h>
#ifdef DEBUG2
#include <cstdio>
#endif
class Interpolation
{
public:

    Vector3 ****pSurfaceInterpolation; ///< surface interpolation polinomials

    unsigned int grade;                ///< the grade aka number of neighbours that constitute one polinomial
    unsigned int noPatch;              ///< number of patches
    unsigned int noHlpElements;        ///< number of refinements

    double h;                          ///< 1/step_size on reference domain

    Interpolation() { grade = 0; noHlpElements = 0; h = 0;};

    /**
     * @brief computes the coefficients of the interpolation polynomial
     *
     * @param[in] U the points of the surface
     * @param[in] gradeIn the grade for the tensor interpolation - aka the grade
     * in one direction
     * @param[in] type a parameter to define the interpolation type, not used
     * here, in this case only the newton interpolation is used
     * @param[in] noHlpElementsIn number of refinements
     * @param[in] noPatchIn number of patches used
     */
    Interpolation(Vector3*** U, int gradeIn, const int type,
                  unsigned int noHlpElementsIn, unsigned int noPatchIn);

    /**
     * @brief computes the interpolation in the 2 point vector
     *
     * @param[in] a the point for which the surface point needs to be computed
     * @param[in] patch the patch to which the point belongs
     */
    Vector3 Chi(Vector2 a, int patch);

    /**
     * @brief computes the derivative of the interpolation in the 2 point vector
     *
     * @param[in] a the point for which the surface point needs to be computed
     * @param[in] patch the patch to which the point belongs
     */
    Vector3 n_Chi(Vector2 a, int patch);

    /**
     * @brief the destructor for the interpolation class, deallocates the
     * interpolation polinomials
     */
    ~Interpolation();
};
#endif

