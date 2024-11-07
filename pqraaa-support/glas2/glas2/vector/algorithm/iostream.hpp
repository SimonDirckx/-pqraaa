//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_algorithm_iostream_hpp
#define glas2_vector_algorithm_iostream_hpp

#include <glas2/vector/concept/vector.hpp>
#include <glas2/concept/is.hpp>
#include <iostream>
#include <type_traits>

namespace glas2 {

  template <typename V>
  typename std::enable_if< is<Vector,V>::value, std::ostream&>::type operator<< ( std::ostream& s, V const& v ) {
    s << "(" << v.size() << ")[" ;
    if (v.size()>0) s << v(0) ;
    for (typename V::size_type i=1; i<v.size(); ++i) { s << "," << v(i) ; }
    s << "]" ;
    return s ;
  }

} // namespace glas2

#endif
