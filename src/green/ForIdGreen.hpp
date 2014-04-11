#ifndef FORIDGREEN_HPP
#define FORIDGREEN_HPP

#include <stdexcept>
#include <string>

#include "Config.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/mpl/begin_end.hpp>
#include <boost/mpl/deref.hpp>
#include <boost/mpl/distance.hpp>
#include <boost/mpl/find.hpp>
#include <boost/mpl/next_prior.hpp>
#include <boost/mpl/vector.hpp>

#include "GreenData.hpp"
#include "IGreensFunction.hpp"

/*! \file ForIdGreen.hpp
 *
 *  Metafunction for fake run-time selection of template arguments.
 *  D. Langr, P.Tvrdik, T. Dytrych and J. P. Draayer
 *  "Fake Run-Time Selection of Template Arguments in C++"
 *  in "Objects, Models, Components, Patterns"
 *  DOI: 10.1007/978-3-642-30561-0_11
 *  http://arxiv.org/abs/1306.5142
 *
 *  The idea is to iterate over a sequence of types, S, beginning
 *  at type B and ending at type E. This heavily relies on Boost MPL.
 *  In the original implementation this is done for an arbitrary
 *  number of template parameters, here we use the 1-dimensional
 *  version of the original code.
 *  Once the template parameter has been resolved, a functor is
 *  applied. The return value of the functor is the Green's function
 *  we want to create with the template parameter resolved.
 *
 *  MPL metafunctions used:
 *  . boost::mpl::distance (http://www.boost.org/doc/libs/1_55_0b1/libs/mpl/doc/refmanual/distance.html); 
 *  . boost::mpl::begin (http://www.boost.org/doc/libs/1_55_0b1/libs/mpl/doc/refmanual/begin.html);
 *  . boost::mpl::end (http://www.boost.org/doc/libs/1_55_0b1/libs/mpl/doc/refmanual/end.html);
 *  . boost::mpl::find (http://www.boost.org/doc/libs/1_55_0b1/libs/mpl/doc/refmanual/find.html);
 *  . boost::mpl::deref (http://www.boost.org/doc/libs/1_55_0b1/libs/mpl/doc/refmanual/deref.html);
 *  . boost::mpl::next (http://www.boost.org/doc/libs/1_55_0b1/libs/mpl/doc/refmanual/next.html);
 *
 */

namespace mpl = boost::mpl;


/*! \struct pos
 *  \brief Returns a zero-based index of a type within a type sequence
 *  \tparam S type sequence
 *  \tparam T the type whose index will be returned
 */
template <typename S, typename T>
struct pos : mpl::distance<
        typename mpl::begin<S>::type,
        typename mpl::find<S, T>::type
        >::type { };

// Primary template
/*! \struct for_id_impl_1
 *  \tparam S type sequence
 *  \tparam B type of the first element in S
 *  \tparam E type of the last element in S
 */
template <
typename S,
         typename B = typename mpl::begin<S>::type,
         typename E = typename mpl::end<S>::type
         >
struct for_id_impl_1 {
    /*! \fn template <typename T> static IGreensFunction * execute(T & f, const greenData & _data, int id)
     *  \brief Iterates over a type sequence either until the position of the actual type matches the
     *         desired id of until the end of the sequence is reached.
     *  \param     f the creational functor to be applied
     *  \param _data the data needed to create the correct Green's function
     *  \param    id the type of derivative
     *  \tparam    T the type of the Green's function
     */
    template <typename T>
    static IGreensFunction * execute(T & f, const greenData & _data, int id) {
        if (pos<S, typename mpl::deref<B>::type>::value == id) {
            return (f.template operator()<typename mpl::deref<B>::type>(_data));
        } else if (1 == mpl::distance<B, E>::value) {
            throw std::invalid_argument("Invalid derivative type (id = " +
                                        boost::lexical_cast<std::string>(id) + ") in for_id metafunction.");
        } else {
            return (for_id_impl_1<S, typename mpl::next<B>::type, E>::execute(f, _data, id));
        }
    }
};

// Partial specialization
// It is never reached at run-time, it is needed to stop the recursive instantiation at compile time
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wextra"
template <typename S, typename E>
struct for_id_impl_1<S, E, E> {
    template <typename T>
    static IGreensFunction * execute(T & f, const greenData & _data, int id) { return NULL; }
};
#pragma GCC diagnostic pop

// Wrapper to the primary template
template <typename S, typename T>
IGreensFunction * for_id(T & f, const greenData & _data, int id)
{
    return (for_id_impl_1<S>::execute(f, _data, id));
}

#endif // FORIDGREEN_HPP
