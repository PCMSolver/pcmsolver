/**
 * @file Topology.hpp
 *
 * @brief computes topological data
 */
#ifndef _TOPOLOGY_HPP_
#define _TOPOLOGY_HPP_

class Vector3;

/// calculates the union K(d,r) = K(d1,r1) union K(d2,r2)
void unify(Vector3 *d, double *r, Vector3 d1, double r1, Vector3 d2, double r2);

#endif
