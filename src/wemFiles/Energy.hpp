/**
 * @file Energy.hpp
 *
 * @brief calculates the potential energy
 */

#ifndef ENERGY_HPP
#define ENERGY_HPP

class GenericAnsatzFunction;

double energy(double *u, GenericAnsatzFunction *af);

double energy_ext(double *u, double *potential, GenericAnsatzFunction *af);

double charge_ext(double * u, double *charge, GenericAnsatzFunction *af);
#endif
