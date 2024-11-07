//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#include <glas2/concept/traits.hpp>

namespace glas2 {

  template <typename X>
  traits<X> info( X& ) {
    return traits<X>() ;
  } //info()

} // namespace glas2
