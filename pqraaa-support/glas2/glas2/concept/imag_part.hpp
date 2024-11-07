//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_concept_imag_part_hpp
#define glas2_concept_imag_part_hpp

#include <type_traits>
#include <glas2/concept/imag.hpp>
#include <glas2/type/imag_proxy.hpp>

namespace glas2 {

  template <typename X>
  struct imag_part {
    typename X::value_type operator() ( X const& x) const { return x.imag() ; }
    //imag_proxy<X> operator() ( X const& x) const { return imag_proxy<X>(x) ; }

    imag_proxy<X> operator() ( X& x) const { return imag_proxy<X>(x) ; }

    typedef imag_proxy<X> result_type ;
  } ;

}

#endif
