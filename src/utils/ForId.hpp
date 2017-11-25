/*
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2017 Roberto Di Remigio, Luca Frediani and collaborators.
 *
 * This file is part of PCMSolver.
 *
 * PCMSolver is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PCMSolver is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to the
 * PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
 */

#pragma once

#include <stdexcept>
#include <string>

#include "Config.hpp"

#include <boost/mpl/begin_end.hpp>
#include <boost/mpl/deref.hpp>
#include <boost/mpl/distance.hpp>
#include <boost/mpl/find.hpp>
#include <boost/mpl/next_prior.hpp>
#include <boost/mpl/vector.hpp>

/*! \file ForId.hpp
 *  Metafunction for fake run-time selection of template arguments.
 *  D. Langr, P.Tvrdik, T. Dytrych and J. P. Draayer
 *  "Fake Run-Time Selection of Template Arguments in C++"
 *  in "Objects, Models, Components, Patterns"
 *  DOI: 10.1007/978-3-642-30561-0_11
 *  http://arxiv.org/abs/1306.5142
 *
 *  The idea is to iterate over sequences of types, S, beginning
 *  at type B and ending at type E. This heavily relies on the Boost
 *  MetaProgramming Library (Boost.MPL)
 *  We provide here an implementation covering the case were up to
 *  three template parameters have to be selected.
 *  According to the original reference, if D is the maximum number
 *  of template parameters to be selected at run-time, then the metaprogramming
 *  technique here described requires coding of (2D + 1) structs:
 *      - one primary template. The primary template has (4D + 1)
 *      template parameters, 4D of which type template.
 *      The primary template exhausts the first type sequence;
 *      - a first partial specialization, exhausting the second
 *      type sequence. This partial specialization has 4D template parameters,
 *      (4D - 1) of which type template;
 *      - a second partial specialization, exhausting the third type
 *      sequence. This partial specialization has (4D - 1) template parameters,
 *      (4D - 2) of which type templates;
 *      - a third partial specialization, never reached at run-time.
 *      This is needed to stop recursive instantiation at compile-time.
 *      It has (4D - 2) template parameters, (4D - 3) of which type template.
 *      - D ApplyFunctor helper structs, needed to wrap the functor.
 *  Once the template parameter has been resolved, a functor is
 *  applied. The return value of the functor is the time-dependent solver
 *  we want to create with the template parameters resolved.
 *
 *  MPL metafunctions used:
 *      - boost::mpl::distance
 *      - boost::mpl::begin
 *      - boost::mpl::end
 *      - boost::mpl::find
 *      - boost::mpl::deref
 *      - boost::mpl::next
 */

namespace mpl = boost::mpl;

/*! \brief Returns a zero-based index of a type within a type sequence
 *  \tparam S type sequence
 *  \tparam T the type whose index will be returned
 */
template <typename S, typename T>
struct position : mpl::distance<typename mpl::begin<S>::type,
                                typename mpl::find<S, T>::type>::type {};

/*! Convenience wrapper to the functor
 *  \tparam D  dimension of the problem, i.e. the maximum number of type sequences
 * allowed
 *  \tparam T1 type selected from S1
 *  \tparam T2 type selected from S2
 *  \tparam T3 type selected from S3
 */
template <int D, typename T1, typename T2, typename T3> struct ApplyFunctor;

/*! \brief Iterates over type sequences either until the position of the actual type
 *matches the
 *         desired id or until the end of the sequences is reached.
 *  \tparam D  dimension of the problem, i.e. the maximum number of type sequences
 *allowed
 *  \tparam S1 type sequence number 1
 *  \tparam S2 type sequence number 2, by default empty
 *  \tparam S3 type sequence number 3, by default empty
 *  \tparam B1 type of the first element in S1
 *  \tparam B2 type of the first element in S2
 *  \tparam B3 type of the first element in S3
 *  \tparam E1 type of the last element in S1
 *  \tparam E2 type of the last element in S2
 *  \tparam E3 type of the last element in S3
 *  \tparam T1 type selected from S1
 *  \tparam T2 type selected from S2
 *  \tparam T3 type selected from S3
 *
 *  This is the primary template. As such, it has (4D + 1) template parameters, 4D
 *  of which type template parameters.
 *
 *  \note The mechanics of the code is also documented, since this is quite an
 *intricate piece of code!
 */
template <int D,
          typename S1,
          typename S2 = mpl::vector<>,
          typename S3 = mpl::vector<>,
          typename B1 = typename mpl::begin<S1>::type,
          typename B2 = typename mpl::begin<S2>::type,
          typename B3 = typename mpl::begin<S3>::type,
          typename E1 = typename mpl::end<S1>::type,
          typename E2 = typename mpl::end<S2>::type,
          typename E3 = typename mpl::end<S3>::type,
          typename T1 = typename mpl::deref<B1>::type,
          typename T2 = typename mpl::deref<B2>::type,
          typename T3 = typename mpl::deref<B3>::type>
struct for_id_impl {
  /*! Executes given functor with selected template arguments
   *  \param[in] f the creational functor to be applied
   *  \param[in] data the data needed to create the correct time-dependent solver
   *  \param[in] id1 the type of derivative
   *  \param[in] id2 index for the second type
   *  \param[in] id3 index for the third type
   *  \tparam ReturnType type of the object returned
   *  \tparam T the type of the functor
   *  \tparam InputType type of the input data
   */
  template <typename ReturnType, typename T, typename InputType>
  static ReturnType * execute(T & f,
                              const InputType & data,
                              int id1,
                              int id2 = 0,
                              int id3 = 0) {
    if (position<S1, typename mpl::deref<B1>::type>::value ==
        id1) {      // Desired type in S1 found
      if (1 == D) { // One-dimensional, we're done!
        return ApplyFunctor<D, typename mpl::deref<B1>::type, T2, T3>::
            template apply<ReturnType>(f, data);
      } else { // Resolve second and third dimensions
        // Call first partial specialization of for_id_impl
        // B1 is not "passed", since the type desired from S1 has been resolved
        // The resolved type is "saved" into T1 and "passed" to the first partial
        // specialization
        return (
            for_id_impl<D,
                        S1,
                        S2,
                        S3,
                        E1,
                        B2,
                        B3,
                        E1,
                        E2,
                        E3,
                        typename mpl::deref<B1>::type,
                        T2>::template execute<ReturnType>(f, data, id1, id2, id3));
      }
    } else if (1 == mpl::distance<B1, E1>::value) { // Desired type NOT found in S1
      throw std::invalid_argument("Invalid derivative type (id1 = " +
                                  pcm::to_string(id1) + ") in for_id metafunction.");
    } else { // First type not resolved, but S1 type sequence not exhausted
      // Call for_id_impl primary template with type of B1 moved to the next type in
      // S1
      return (for_id_impl<D,
                          S1,
                          S2,
                          S3,
                          typename mpl::next<B1>::type,
                          B2,
                          B3,
                          E1,
                          E2,
                          E3,
                          T1,
                          T2,
                          T3>::template execute<ReturnType>(f, data, id1, id2, id3));
    }
  }
};

/*! \brief Iterates over type sequences either until the position of the actual type
 *matches the
 *         desired id or until the end of the sequences is reached.
 *  \tparam D  dimension of the problem, i.e. the maximum number of type sequences
 *allowed
 *  \tparam S1 type sequence number 1
 *  \tparam S2 type sequence number 2
 *  \tparam S3 type sequence number 3
 *  \tparam B2 type of the first element in S2
 *  \tparam B3 type of the first element in S3
 *  \tparam E1 type of the last element in S1
 *  \tparam E2 type of the last element in S2
 *  \tparam E3 type of the last element in S3
 *  \tparam T1 type selected from S1
 *  \tparam T2 type selected from S2
 *  \tparam T3 type selected from S3
 *
 *  This is the fist partial specialization. As such, it has 4D template parameters,
 *(4D - 1)
 *  of which type template parameters.
 *
 *  \note The mechanics of the code is also documented, since this is quite an
 *intricate piece of code!
 */
template <int D,
          typename S1,
          typename S2,
          typename S3,
          typename B2,
          typename B3,
          typename E1,
          typename E2,
          typename E3,
          typename T1,
          typename T2,
          typename T3>
struct for_id_impl<D, S1, S2, S3, E1, B2, B3, E1, E2, E3, T1, T2, T3> {
  /*! Executes given functor with selected template arguments
   *  \param[in] f the creational functor to be applied
   *  \param[in] data the data needed to create the correct time-dependent solver
   *  \param[in] id1 the type of derivative
   *  \param[in] id2 index for the second type
   *  \param[in] id3 index for the third type
   *  \tparam ReturnType type of the object returned
   *  \tparam T the type of the time-dependent solver
   *  \tparam InputType type of the input data
   */
  template <typename ReturnType, typename T, typename InputType>
  static ReturnType * execute(T & f,
                              const InputType & data,
                              int id1,
                              int id2 = 0,
                              int id3 = 0) {
    if (position<S2, typename mpl::deref<B2>::type>::value ==
        id2) {      // Desired type in S2 found
      if (2 == D) { // Two-dimensional, we're done!
        return ApplyFunctor<D, T1, typename mpl::deref<B2>::type, T3>::
            template apply<ReturnType>(f, data);
      } else { // Resolve third dimension
        // Call second partial specialization of for_id_impl
        // B2 is not "passed", since the type desired from S2 has been resolved
        // The resolved type is "saved" into T2 and "passed" to the first partial
        // specialization
        return (for_id_impl<D,
                            S1,
                            S2,
                            S3,
                            E1,
                            E2,
                            B3,
                            E1,
                            E2,
                            E3,
                            typename mpl::deref<B2>::type>::
                    template execute<ReturnType>(f, data, id1, id2, id3));
      }
    } else if (1 == mpl::distance<B2, E2>::value) { // Desired type NOT found in S2
      throw std::invalid_argument("Invalid integrator policy (id2 = " +
                                  pcm::to_string(id2) + ") in for_id metafunction.");
    } else { // Second type not resolved, but S2 type sequence not exhausted
      // Call for_id_impl first partial specialization with type of B2 moved to the
      // next type in S2
      return (for_id_impl<D,
                          S1,
                          S2,
                          S3,
                          E1,
                          typename mpl::next<B2>::type,
                          B3,
                          E1,
                          E2,
                          E3,
                          T1,
                          T2,
                          T3>::template execute<ReturnType>(f, data, id1, id2, id3));
    }
  }
};

/*! \brief Iterates over type sequences either until the position of the actual type
 *matches the
 *         desired id or until the end of the sequences is reached.
 *  \tparam D  dimension of the problem, i.e. the maximum number of type sequences
 *allowed
 *  \tparam S1 type sequence number 1
 *  \tparam S2 type sequence number 2
 *  \tparam S3 type sequence number 3
 *  \tparam B3 type of the first element in S3
 *  \tparam E1 type of the last element in S1
 *  \tparam E2 type of the last element in S2
 *  \tparam E3 type of the last element in S3
 *  \tparam T1 type selected from S1
 *  \tparam T2 type selected from S2
 *  \tparam T3 type selected from S3
 *
 *  This is the second partial specialization. As such, it has (4D - 1) template
 *parameters, (4D - 2)
 *  of which type template parameters.
 *
 *  \note The mechanics of the code is also documented, since this is quite an
 *intricate piece of code!
 */
template <int D,
          typename S1,
          typename S2,
          typename S3,
          typename B3,
          typename E1,
          typename E2,
          typename E3,
          typename T1,
          typename T2,
          typename T3>
struct for_id_impl<D, S1, S2, S3, E1, E2, B3, E1, E2, E3, T1, T2, T3> {
  /*! Executes given functor with selected template arguments
   *  \param[in] f the creational functor to be applied
   *  \param[in] data the data needed to create the correct time-dependent solver
   *  \param[in] id1 the type of derivative
   *  \param[in] id2 index for the second type
   *  \param[in] id3 index for the third type
   *  \tparam ReturnType type of the object returned
   *  \tparam T the type of the time-dependent solver
   *  \tparam InputType type of the input data
   */
  template <typename ReturnType, typename T, typename InputType>
  static ReturnType * execute(T & f,
                              const InputType & data,
                              int id1,
                              int id2 = 0,
                              int id3 = 0) {
    if (position<S3, typename mpl::deref<B3>::type>::value ==
        id3) { // Desired type in S3 found, we're done!
      return ApplyFunctor<D, T1, T2, typename mpl::deref<B3>::type>::template apply<
          ReturnType>(f, data);
    } else if (1 == mpl::distance<B3, E3>::value) { // Desired type NOT found in S3
      throw std::invalid_argument("Invalid permittivity profile (id3 = " +
                                  pcm::to_string(id3) + ") in for_id metafunction.");
    } else { // Third type not resolved, but S3 type sequence not exhausted
      // Call for_id_impl second partial specialization with type of B3 moved to the
      // next type in S3
      return (for_id_impl<D,
                          S1,
                          S2,
                          S3,
                          typename mpl::next<B3>::type,
                          E1,
                          E2,
                          E3,
                          T1,
                          T2,
                          T3>::template execute<ReturnType>(f, data, id1, id2, id3));
    }
  }
};

/*! \brief Iterates over type sequences either until the position of the actual type
 *matches the
 *         desired id or until the end of the sequences is reached.
 *  \tparam D  dimension of the problem, i.e. the maximum number of type sequences
 *allowed
 *  \tparam S1 type sequence number 1
 *  \tparam S2 type sequence number 2
 *  \tparam S3 type sequence number 3
 *  \tparam E1 type of the last element in S1
 *  \tparam E2 type of the last element in S2
 *  \tparam E3 type of the last element in S3
 *  \tparam T1 type selected from S1
 *  \tparam T2 type selected from S2
 *  \tparam T3 type selected from S3
 *
 *  This is the third partial specialization. As such, it has (4D - 2) template
 *parameters, (4D - 3)
 *  of which type template parameters.
 *  It is never reached at run-time, it is needed to stop the recursive instantiation
 *at compile-time.
 */
template <int D,
          typename S1,
          typename S2,
          typename S3,
          typename E1,
          typename E2,
          typename E3,
          typename T1,
          typename T2,
          typename T3>
struct for_id_impl<D, S1, S2, S3, E1, E2, E3, E1, E2, E3, T1, T2, T3> {
  /*! Executes given functor with selected template arguments
   *  \tparam ReturnType type of the object returned
   *  \tparam T the type of the time-dependent solver
   *  \tparam InputType type of the input data
   */
  template <typename ReturnType, typename T, typename InputType>
  static ReturnType * execute(T & /* f */,
                              const InputType & /* data */,
                              int /* id1 */,
                              int /* id2 */ = 0,
                              int /* id3 */ = 0) {
    return __nullptr;
  }
};

/*! @{ Wrappers to the functor */
/*! Wrapper to the three-dimensional case
 *  \tparam T1 type selected from S1
 *  \tparam T2 type selected from S2
 *  \tparam T3 type selected from S3
 */
template <typename T1, typename T2, typename T3> struct ApplyFunctor<3, T1, T2, T3> {
  /*! Executes given functor with selected template arguments
   *  \tparam ReturnType type of the object returned
   *  \tparam T the type of the time-dependent solver
   *  \tparam InputType type of the input data
   */
  template <typename ReturnType, typename T, typename InputType>
  static ReturnType * apply(T & f, const InputType & data) {
    return (f.template operator()<T1, T2, T3>(data));
  }
};

/*! Wrapper to the two-dimensional case
 *  \tparam T1 type selected from S1
 *  \tparam T2 type selected from S2
 *  \tparam T3 type selected from S3
 */
template <typename T1, typename T2, typename T3> struct ApplyFunctor<2, T1, T2, T3> {
  /*! Executes given functor with selected template arguments
   *  \tparam ReturnType type of the object returned
   *  \tparam T the type of the time-dependent solver
   *  \tparam InputType type of the input data
   */
  template <typename ReturnType, typename T, typename InputType>
  static ReturnType * apply(T & f, const InputType & data) {
    return (f.template operator()<T1, T2>(data));
  }
};

/*! Wrapper to the one-dimensional case
 *  \tparam T1 type selected from S1
 *  \tparam T2 type selected from S2
 *  \tparam T3 type selected from S3
 */
template <typename T1, typename T2, typename T3> struct ApplyFunctor<1, T1, T2, T3> {
  /*! Executes given functor with selected template arguments
   *  \tparam ReturnType type of the object returned
   *  \tparam T the type of the time-dependent solver
   *  \tparam InputType type of the input data
   */
  template <typename ReturnType, typename T, typename InputType>
  static ReturnType * apply(T & f, const InputType & data) {
    return (f.template operator()<T1>(data));
  }
};
/*! @}*/

/*! @{ Wrappers to the primary template for the metafunction */
/*! Wrapper for the three-dimensional case.
 *  \tparam S1 type sequence 1
 *  \tparam S2 type sequence 2
 *  \tparam S3 type sequence 3
 *  \tparam ReturnType type of the object returned
 *  \tparam T  type of the functor
 *  \tparam InputType type of the input data
 *  \param[in,out] f functor
 *  \param[in] data input arguments for functor
 *  \param[in] id1  index of the first template argument, selected in S1
 *  \param[in] id2  index of the second template argument, selected in S2
 *  \param[in] id3  index of the third template argument, selected in S3
 */
template <typename S1,
          typename S2,
          typename S3,
          typename ReturnType,
          typename T,
          typename InputType>
ReturnType * for_id(T & f, const InputType & data, int id1, int id2, int id3) {
  return for_id_impl<3, S1, S2, S3>::template execute<ReturnType>(
      f, data, id1, id2, id3);
}

/*! Wrapper for the two-dimensional case.
 *  \tparam S1 type sequence 1
 *  \tparam S2 type sequence 2
 *  \tparam ReturnType type of the object returned
 *  \tparam T  type of the functor
 *  \tparam InputType type of the input data
 *  \param[in,out] f functor
 *  \param[in] data input arguments for functor
 *  \param[in] id1  index of the first template argument, selected in S1
 *  \param[in] id2  index of the second template argument, selected in S2
 */
template <typename S1,
          typename S2,
          typename ReturnType,
          typename T,
          typename InputType>
ReturnType * for_id(T & f, const InputType & data, int id1, int id2) {
  return for_id_impl<2, S1, S2>::template execute<ReturnType>(f, data, id1, id2);
}

/*! Wrapper for the one-dimensional case.
 *  \tparam S1 type sequence 1
 *  \tparam ReturnType type of the object returned
 *  \tparam T  type of the functor
 *  \tparam InputType type of the input data
 *  \param[in,out] f functor
 *  \param[in] data input arguments for functor
 *  \param[in] id1  index of the first template argument, selected in S1
 */
template <typename S1, typename ReturnType, typename T, typename InputType>
ReturnType * for_id(T & f, const InputType & data, int id1) {
  return for_id_impl<1, S1>::template execute<ReturnType>(f, data, id1);
}
/*! @}*/
