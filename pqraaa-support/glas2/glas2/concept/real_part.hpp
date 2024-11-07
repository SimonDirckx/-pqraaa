//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_concept_real_part_hpp
#define glas2_concept_real_part_hpp

#include <type_traits>
#include <glas2/concept/real.hpp>
#include <glas2/type/real_proxy.hpp>

namespace glas2 {

  struct real_part {
    template <typename X>
    typename X::value_type operator() ( X const& x) const { return x.real() ; }
    //real_proxy<X> operator() ( X const& x) const { return real_proxy<X>(x) ; }

    template <typename X>
    real_proxy<X> operator() ( X& x) const { return real_proxy<X>(x) ; }
  } ;

}

#endif
