#ifndef FORID_HPP
#define FORID_HPP

#include <stdexcept>
#include <string>

#include <boost/lexical_cast.hpp>
#include <boost/mpl/begin_end.hpp>
#include <boost/mpl/deref.hpp>
#include <boost/mpl/distance.hpp>
#include <boost/mpl/find.hpp>
#include <boost/mpl/next_prior.hpp>
#include <boost/mpl/vector.hpp>

/*! \file for_id.hpp
 *
 *  Metafunction for fake run-time selection of template arguments.
 *  D. Langr, P. Tvrdik, T. Dytrych and J. P. Draayer
 *  "Fake Run-Time Selection of Template Arguments in C++"
 *  in "Objects, Models, Components, Patterns"
 *  DOI: 10.1007/978-3-642-30561-0_11
 *  http://arxiv.org/abs/1306.5142
 *
 *  The idea is to iterate over a sequence of types, S, beginning
 *  at type B and ending at type E. This heavily relies on Boost MPL.
 *  In the original implementation this is done for an arbitrary
 *  number of template parameters.
 *  Once the template parameter has been resolved, a functor is
 *  applied.
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

template <typename S, typename T>
struct pos : mpl::distance<
        typename mpl::begin<S>::type,
        typename mpl::find<S, T>::type
        >::type { };

template <int D, typename T1, typename T2>
struct executor;

template <typename T1, typename T2>
struct executor<1, T1, T2> {
    template <typename T>
    static void execute(T & f) {
        f.template operator()<T1>();
    }
};

template <typename T1, typename T2>
struct executor<2, T1, T2> {
    template <typename T>
    static void execute(T & f) {
        f.template operator()<T1, T2>();
    }
};


// Primary template
template <
int D,
    typename S1,
    typename S2 = mpl::vector<>,
    typename B1 = typename mpl::begin<S1>::type,
    typename B2 = typename mpl::begin<S2>::type,
    typename E1 = typename mpl::end<S1>::type,
    typename E2 = typename mpl::end<S2>::type,
    typename T1 = typename mpl::deref<B1>::type,
    typename T2 = typename mpl::deref<B2>::type
    >
struct for_id_impl {
    template <typename T>
    static void execute(T & f, int id1, int id2 = 0) {
        if (pos<S1, typename mpl::deref<B1>::type>::value == id1)
            if (1 == D)
                executor<D, typename mpl::deref<B1>::type, T2>::execute(f);
            else
                for_id_impl<
                D, S1, S2, E1, B2, E1, E2, typename mpl::deref<B1>::type
                >::execute(f, id1, id2);
        else if (1 == mpl::distance<B1, E1>::value)
            throw std::invalid_argument("Invalid type (id = " +
                                        boost::lexical_cast<std::string>(id) + ") in for_id metafunction.");
        else
            for_id_impl<
            D, S1, S2, typename mpl::next<B1>::type, B2, E1, E2, T1
            >::execute(f, id1);
    }
};

// Partial specialization #1
template <
int D,
    typename S1,
    typename S2,
    typename B2,
    typename E1,
    typename E2,
    typename T1,
    typename T2
    >
struct for_id_impl<D, S1, S2, E1, B2, E1, E2, T1, T2> {
    template <typename T>
    static void execute(T & f, int id1, int id2 = 0) {
        if (pos<S2, typename mpl::deref<B2>::type>::value == id2)
            executor<D, T1, typename mpl::deref<B2>::type>::execute(f);
        else if (1 == mpl::distance<B2, E2>::value)
            throw std::invalid_argument("Invalid type (id = " +
                                        boost::lexical_cast<std::string>(id) + ") in for_id metafunction.");
        else
            for_id_impl<
            D, S1, S2, E1, typename mpl::next<B2>::type, E1, E2, T1
            >::execute(f, id1, id2);
    }
};

// Partial specialization #2
template <
int D,
    typename S1,
    typename S2,
    typename E1,
    typename E2,
    typename T1,
    typename T2
    >
struct for_id_impl<D, S1, S2, E1, E2, E1, E2, T1, T2> {
    template <typename T>
    static void execute(T & f, int id1, int id2 = 0) { }
};

// One-dimensional case
template <typename S1, typename T>
void for_id(T & f, int id1)
{
    for_id_impl<1, S1>::execute(f, id1);
}

// Two-dimensional case
template <typename S1, typename S2, typename T>
void for_id(T & f, int id1, int id2)
{
    for_id_impl<2, S1, S2>::execute(f, id1, id2);
}

#endif // FORID_HPP
