//Zifan Liu zifan.liu@cs.kuleuven.be

#ifndef glas2_bsp_vector_algorithm_inner_prod_hpp
#define glas2_bsp_vector_algorithm_inner_prod_hpp

#include <glas2/bsp/backend/current_backend.hpp>
#include <glas2/bsp/backend/default/vector/inner_prod.hpp>
#include <glas2/bsp/vector/concept/bsp_vector.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>


namespace glas2 { namespace bsp {

  template <typename X, typename Y>
  typename std::enable_if< is<BSPVector,X>::value && is<BSPVector,Y>::value
                         , decltype( typename X::value_type() * typename Y::value_type() )
                         >::type inner_prod( X const& x, Y const& y ) {
	   

     return glas2::bsp::inner_prod( bsp::current_backend(), x, y ) ;


  }
} } // namespace glas2::bsp

#endif
