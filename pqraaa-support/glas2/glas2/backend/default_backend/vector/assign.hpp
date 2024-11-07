//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_backend_default_backend_vector_assign_hpp
#define glas2_backend_default_backend_vector_assign_hpp

#include <glas2/backend/default_backend/default_backend.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename To, typename From>
  typename std::enable_if< is<DenseVector,To>::value && is<DenseVector,From>::value, To& >::type assign( default_backend, To& to, From const& from ) {
    assert( to.size() == from.size() ) ;
    for (typename To::size_type i=0; i<to.size(); ++i) { to(i) = from(i) ; }
    return to ;
  }

} // namespace glas2

#endif
