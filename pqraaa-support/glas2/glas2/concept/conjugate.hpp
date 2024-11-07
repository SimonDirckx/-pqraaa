//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_concept_conjugate_hpp
#define glas2_concept_conjugate_hpp

#include <type_traits>
#include <glas2/concept/conj.hpp>

namespace glas2 {

  struct conjugate {
    template <typename X>
    decltype(glas2::conj(X())) operator() ( X const& x) const { return glas2::conj(x) ; }
  } ;

}

#endif
