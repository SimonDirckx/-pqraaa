//Zifan Liu zifan.liu@cs.kuleuven.be

#ifndef glas2_bsp_vector_algorithm_norm_hpp
#define glas2_bsp_vector_algorithm_norm_hpp

#include <glas2/bsp/backend/current_backend.hpp>
#include <glas2/bsp/backend/default/vector/norm_2.hpp>
#include <glas2/bsp/vector/concept/bsp_vector.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cmath>


namespace glas2 { namespace bsp {

  template <typename X>
  typename std::enable_if< is<BSPVector,X>::value
                         , decltype( std::abs( typename X::value_type() ) )
                         >::type norm_2( X const& x ) {
	   

     return glas2::bsp::norm_2( bsp::current_backend(), x ) ;


  }
} } // namespace glas2::bsp

#endif
