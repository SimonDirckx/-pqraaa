//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_concept_real_hpp
#define glas2_concept_real_hpp

#include <type_traits>

#ifndef GLAS_COMPLEX
#define GLAS_COMPLEX
#endif

#ifdef GLAS_COMPLEX
#include <complex>
#include <glas2/type/real_proxy.hpp>
#endif

namespace glas2 {

  template <typename X>
  typename std::enable_if< std::is_arithmetic<X>::value, X >::type real( X const& x ) { return x ; }

#ifdef GLAS_COMPLEX
  template <typename X>
  typename std::enable_if< std::is_arithmetic<X>::value, X >::type real( std::complex<X> const& x ) { return x.real() ; }

  template <typename X>
  typename std::enable_if< std::is_arithmetic<X>::value, real_proxy< std::complex<X> > >::type real( std::complex<X>& x ) { return real_proxy< std::complex<X> >(x) ; }
#endif

}

#endif
