#include <stdexcept>

#include <boost/mpl/begin_end.hpp>
#include <boost/mpl/deref.hpp>
#include <boost/mpl/distance.hpp>
#include <boost/mpl/find.hpp>
#include <boost/mpl/next_prior.hpp>
#include <boost/mpl/vector.hpp>

using namespace boost::mpl;

template <typename S, typename T>
struct pos : distance<
	     typename begin<S>::type,
	     typename find<S, T>::type
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
   typename S2 = vector<>,
   typename B1 = typename begin<S1>::type,
   typename B2 = typename begin<S2>::type,
   typename E1 = typename end<S1>::type,
   typename E2 = typename end<S2>::type,
   typename T1 = typename deref<B1>::type,
   typename T2 = typename deref<B2>::type
>
struct for_id_impl {
	template <typename T>
	static void execute(T & f, int id1, int id2 = 0) {
          if (pos<S1, typename deref<B1>::type>::value == id1)
            if (1 == D)
              executor<D, typename deref<B1>::type, T2>::execute(f);
	    else
              for_id_impl<
		 D, S1, S2, E1, B2, E1, E2, typename deref<B1>::type
	      >::execute(f, id1, id2);
	  else if (1 == distance<B1, E1>::value)
	    throw std::invalid_argument("");
	  else
            for_id_impl<
	       D, S1, S2, typename next<B1>::type, B2, E1, E2, T1
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
      if (pos<S2, typename deref<B2>::type>::value == id2)
	 executor<D, T1, typename deref<B2>::type>::execute(f);
      else if (1 == distance<B2, E2>::value)
	 throw std::invalid_argument("");
      else
	 for_id_impl<
            D, S1, S2, E1, typename next<B2>::type, E1, E2, T1
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
void for_id(T & f, int id1) {
   for_id_impl<1, S1>::execute(f, id1);	
}

// Two-dimensional case
template <typename S1, typename S2, typename T> 
void for_id(T & f, int id1, int id2) {
   for_id_impl<2, S1, S2>::execute(f, id1, id2);	
}
